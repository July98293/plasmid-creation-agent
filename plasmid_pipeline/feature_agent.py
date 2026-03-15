"""
FeatureAgent — resolves biological part sequences from online databases.

Source priority (applied per feature):
  1. NCBI GenBank (nuccore) — primary.
  2. Addgene API              — fallback when NCBI returns nothing valid.

For each field received from the intent agent (promoter, polyA, terminator,
selection_marker, origin_of_replication, tags) the agent:
  1. Cleans the name (strips LLM annotations like "(suggested)").
  2. Builds NCBI search queries dynamically from the name — no static alias map.
  3. For each NCBI result: tries annotated sub-feature extraction first, then
     falls back to the full ORIGIN when the record IS the element.
  4. On NCBI failure: queries Addgene using the appropriate catalog filter
     (promoters=, resistance_marker=, tags=, name=), downloads the first
     available GenBank file, and extracts the annotation with
     PlasmidFeatureExtractor.
  5. Validates sequence length against category-level heuristics.
  6. Returns a ResolvedFeature with name, sequence, source, length, validated.

Kozak (GCCACC) is the only hardcoded sequence — it is a universal 6 bp
consensus with no meaningful standalone database record.
"""
from __future__ import annotations

import json
import os
import re
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Tuple

from logging_utils import get_conversation_logger

from .addgene_client import AddgeneClient
from .models import FeatureInput, FeatureOutput, ResolvedFeature, WarningMessage
from .plasmid_feature_extractor import PlasmidFeatureExtractor


# ---------------------------------------------------------------------------
# Fixed constant
# ---------------------------------------------------------------------------

KOZAK_SEQ = "GCCACC"

# ---------------------------------------------------------------------------
# Category configs — class-level heuristics, not feature-name lookups.
# ---------------------------------------------------------------------------

_CATEGORY: Dict[str, Dict[str, Any]] = {
    "promoter": {
        "gb_types":   {"promoter", "regulatory", "enhancer"},
        "min_bp":     100,
        "max_bp":     3000,
        "output_type": "promoter",
    },
    "polya": {
        "gb_types":   {"polyA_signal", "regulatory"},
        "min_bp":     100,
        "max_bp":     600,
        "output_type": "terminator",
    },
    "terminator": {
        "gb_types":   {"terminator", "regulatory"},
        "min_bp":     50,
        "max_bp":     600,
        "output_type": "terminator",
    },
    "selection_marker": {
        "gb_types":   {"CDS", "gene"},
        "min_bp":     400,
        "max_bp":     2000,
        "output_type": "selection_marker",
    },
    "ori": {
        "gb_types":   {"rep_origin", "misc_feature"},
        "min_bp":     200,
        "max_bp":     2500,
        "output_type": "origin_of_replication",
    },
    "tag": {
        "gb_types":   {"CDS", "misc_feature"},
        "min_bp":     18,
        "max_bp":     3000,
        "output_type": "tag",
    },
}

_MAMMALIAN_HOST_TOKENS = ["hek", "293", "cho", "cos", "hela", "mamm", "human", "mouse"]

# ---------------------------------------------------------------------------
# Known fixed NCBI accessions for well-characterised parts.
# These bypass the dynamic title search entirely — the accession is definitive.
# key: substring(s) to match against the cleaned, lowercased feature name.
# ---------------------------------------------------------------------------

_KNOWN_ACCESSIONS: Dict[str, Dict[str, Any]] = {
    # TEM-1 beta-lactamase — the canonical AmpR gene.
    # V00613: Sutcliffe 1979, complete bla CDS, 861 bp.
    "ampicillin": {
        "accession":   "V00613",
        "label":       "TEM-1 beta-lactamase (AmpR)",
        "min_bp":      800,
        "max_bp":      900,
    },
    "ampr": {
        "accession":   "V00613",
        "label":       "TEM-1 beta-lactamase (AmpR)",
        "min_bp":      800,
        "max_bp":      900,
    },
    "tem-1": {
        "accession":   "V00613",
        "label":       "TEM-1 beta-lactamase (AmpR)",
        "min_bp":      800,
        "max_bp":      900,
    },
}


class FeatureAgent:
    NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    NCBI_EFETCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    def __init__(
        self,
        *,
        ncbi_api_key: Optional[str] = None,
        ncbi_email: Optional[str] = None,
        user_agent: str = "plasmid-pipeline/0.2",
    ) -> None:
        self._logger = get_conversation_logger()
        self._ncbi_api_key = ncbi_api_key or os.getenv("NCBI_API_KEY", "")
        self._ncbi_email   = ncbi_email   or os.getenv("NCBI_EMAIL", "")
        self._user_agent   = user_agent
        self._extractor    = PlasmidFeatureExtractor()
        self._addgene      = AddgeneClient(user_agent=user_agent)

    # ------------------------------------------------------------------
    # HTTP helpers (mirrors BackboneAgent pattern — synchronous urllib)
    # ------------------------------------------------------------------

    def _http_get_json(self, url: str, params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        full_url = f"{url}?{urllib.parse.urlencode(params)}"
        print(f"[FEATURE HTTP] GET {full_url}", flush=True)
        req = urllib.request.Request(
            full_url,
            headers={"User-Agent": self._user_agent, "Accept": "application/json"},
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                print("[FEATURE HTTP] ✓ JSON", flush=True)
                return data
        except Exception as e:
            print(f"[FEATURE HTTP] ✗ {e}", flush=True)
            return None

    def _http_get_text(self, url: str, params: Dict[str, Any]) -> Optional[str]:
        full_url = f"{url}?{urllib.parse.urlencode(params)}"
        print(f"[FEATURE HTTP] GET {full_url}", flush=True)
        req = urllib.request.Request(
            full_url,
            headers={"User-Agent": self._user_agent, "Accept": "text/plain"},
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                text = resp.read().decode("utf-8")
                print(f"[FEATURE HTTP] ✓ TEXT ({len(text)} bytes)", flush=True)
                return text
        except Exception as e:
            print(f"[FEATURE HTTP] ✗ {e}", flush=True)
            return None

    # ------------------------------------------------------------------
    # NCBI helpers
    # ------------------------------------------------------------------

    def _ncbi_params(self) -> Dict[str, Any]:
        p: Dict[str, Any] = {}
        if self._ncbi_api_key:
            p["api_key"] = self._ncbi_api_key
        if self._ncbi_email:
            p["email"] = self._ncbi_email
        return p

    def _ncbi_search(self, term: str, retmax: int = 5) -> List[str]:
        data = self._http_get_json(self.NCBI_ESEARCH_URL, {
            "db":      "nuccore",
            "term":    term,
            "retmode": "json",
            "retmax":  retmax,
            **self._ncbi_params(),
        })
        return ((data or {}).get("esearchresult") or {}).get("idlist") or []

    def _ncbi_fetch_genbank(self, uid: str) -> Optional[str]:
        return self._http_get_text(self.NCBI_EFETCH_URL, {
            "db":      "nuccore",
            "id":      uid,
            "rettype": "gbwithparts",
            "retmode": "text",
            **self._ncbi_params(),
        })

    def _parse_origin(self, gb_text: str) -> str:
        """Extract raw DNA from the GenBank ORIGIN section."""
        bases: List[str] = []
        in_origin = False
        for line in gb_text.splitlines():
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin:
                if line.startswith("//"):
                    break
                bases.append(re.sub(r"[^acgtACGT]", "", line))
        return "".join(bases).upper()

    # ------------------------------------------------------------------
    # Name cleaning
    # ------------------------------------------------------------------

    def _clean(self, name: str) -> str:
        """Strip LLM-appended annotations like '(suggested)', '(bacterial)', etc."""
        cleaned = re.sub(r"\s*\([^)]*\)\s*$", "", name, flags=re.IGNORECASE)
        return cleaned.strip()

    # ------------------------------------------------------------------
    # NCBI search term construction (no static alias map — name-driven)
    # ------------------------------------------------------------------

    def _ncbi_search_terms(self, clean_name: str, category: str) -> List[str]:
        """Return ordered NCBI query strings from most specific to broadest."""
        n = clean_name
        if category == "promoter":
            return [
                f'"{n} promoter"[Title]',
                f'"{n}"[All Fields] AND "promoter"[Feature Key]',
            ]
        if category == "polya":
            return [
                f'"{n}"[Title]',
                f'"{n}"[All Fields] AND "polyadenylation signal"[Feature Key]',
                f'"{n}"[All Fields] AND polyA[All Fields]',
            ]
        if category == "terminator":
            return [
                f'"{n} terminator"[Title]',
                f'"{n}"[All Fields] AND "terminator"[Feature Key]',
            ]
        if category == "selection_marker":
            return [
                f'"{n} resistance"[Title] AND "complete"[Title]',
                f'"{n}"[All Fields] AND "antibiotic resistance"[All Fields]',
                f'"{n}"[All Fields] AND "resistance gene"[All Fields]',
            ]
        if category == "ori":
            return [
                f'"{n} origin of replication"[Title]',
                f'"{n}"[All Fields] AND "origin of replication"[All Fields]',
                f'"{n} ori"[All Fields]',
            ]
        if category == "tag":
            return [
                f'"{n}"[Title] AND "coding sequence"[Title]',
                f'"{n}"[All Fields] AND "coding sequence"[All Fields]',
            ]
        return [f'"{n}"[Title]']

    # ------------------------------------------------------------------
    # Addgene search — maps category to the right catalog filter
    # ------------------------------------------------------------------

    def _addgene_search_kwargs(self, clean_name: str, category: str) -> List[Dict[str, Any]]:
        """
        Return a list of Addgene search kwarg dicts to try in order.
        Each dict is passed directly to AddgeneClient.search_plasmids(**kwargs).
        """
        n = clean_name
        if category == "promoter":
            return [{"promoters": n, "page_size": 5}]
        if category in ("polya", "terminator"):
            # No dedicated filter — search by name, hoping it appears in plasmid names
            return [{"name": n, "page_size": 5}]
        if category == "selection_marker":
            return [
                {"resistance_marker": n, "page_size": 5},
                {"name": n, "page_size": 5},
            ]
        if category == "ori":
            return [{"name": n, "page_size": 5}]
        if category == "tag":
            return [{"tags": n, "page_size": 5}]
        return [{"name": n, "page_size": 5}]

    # ------------------------------------------------------------------
    # Sequence extraction from a GenBank text blob
    # ------------------------------------------------------------------

    def _extract_seq(
        self,
        gb_text: str,
        clean_name: str,
        cat: Dict[str, Any],
    ) -> Optional[str]:
        """
        Try annotated sub-feature first (PlasmidFeatureExtractor), then fall
        back to the full ORIGIN when the record is a standalone element.
        """
        extracted = self._extractor.extract_by_keywords(
            genbank_text=gb_text,
            desired_feature_types=cat["gb_types"],
            keywords=[clean_name] + clean_name.lower().split(),
        )
        return extracted.sequence if extracted else self._parse_origin(gb_text)

    def _validate_length(
        self,
        seq: str,
        clean_name: str,
        uid_label: str,
        cat: Dict[str, Any],
        warnings: List[WarningMessage],
    ) -> bool:
        length = len(seq)
        if cat["min_bp"] <= length <= cat["max_bp"]:
            return True
        warnings.append(WarningMessage(
            code="FEATURE_LENGTH_MISMATCH",
            message=(
                f"'{clean_name}' from {uid_label}: {length} bp is outside "
                f"expected {cat['min_bp']}–{cat['max_bp']} bp. Skipping."
            ),
        ))
        return False

    # ------------------------------------------------------------------
    # Source 1: NCBI GenBank
    # ------------------------------------------------------------------

    def _fetch_by_known_accession(
        self,
        clean_name: str,
        cat: Dict[str, Any],
    ) -> Optional[ResolvedFeature]:
        """
        For parts with a definitive, fixed NCBI accession: fetch it directly,
        bypassing any title search.  Raises ValueError if the fetch fails or the
        returned sequence is outside the expected length range.
        """
        name_lower = clean_name.lower()
        info = next((v for k, v in _KNOWN_ACCESSIONS.items() if k in name_lower), None)
        if info is None:
            return None

        accession = info["accession"]
        print(
            f"[FEATURE NCBI] '{clean_name}' → fixed accession {accession} ({info['label']})",
            flush=True,
        )

        gb_text = self._ncbi_fetch_genbank(accession)
        if not gb_text:
            raise ValueError(
                f"Failed to fetch NCBI accession {accession} for '{clean_name}' ({info['label']})."
            )

        seq = self._extract_seq(gb_text, info["label"], cat)
        if not seq:
            raise ValueError(
                f"No sequence found in NCBI accession {accession} for '{clean_name}'."
            )

        length = len(seq)
        if not (info["min_bp"] <= length <= info["max_bp"]):
            raise ValueError(
                f"NCBI accession {accession} for '{clean_name}' returned {length} bp; "
                f"expected {info['min_bp']}–{info['max_bp']} bp.  "
                f"The accession may have changed — please verify {accession}."
            )

        print(
            f"[FEATURE NCBI] ✓ '{clean_name}' accession={accession} ({length} bp)",
            flush=True,
        )
        return ResolvedFeature(
            name=info["label"],
            type=cat["output_type"],
            sequence=seq,
            source=f"NCBI:{accession}",
            length=length,
            validated=True,
            position_hint=None,
        )

    def _try_ncbi(
        self,
        clean_name: str,
        category: str,
        cat: Dict[str, Any],
        warnings: List[WarningMessage],
    ) -> Optional[ResolvedFeature]:
        for term in self._ncbi_search_terms(clean_name, category):
            ids = self._ncbi_search(term)
            print(f"[FEATURE NCBI] '{clean_name}' query='{term}' → {ids}", flush=True)
            if not ids:
                continue

            for uid in ids:
                gb_text = self._ncbi_fetch_genbank(uid)
                if not gb_text:
                    continue

                seq = self._extract_seq(gb_text, clean_name, cat)
                if not seq:
                    continue

                if not self._validate_length(seq, clean_name, f"NCBI:{uid}", cat, warnings):
                    continue

                print(f"[FEATURE NCBI] ✓ '{clean_name}' uid={uid} ({len(seq)} bp)", flush=True)
                return ResolvedFeature(
                    name=clean_name,
                    type=cat["output_type"],
                    sequence=seq,
                    source=f"NCBI:{uid}",
                    length=len(seq),
                    validated=True,
                    position_hint=None,
                )

        return None

    # ------------------------------------------------------------------
    # Source 2: Addgene
    # ------------------------------------------------------------------

    def _try_addgene(
        self,
        clean_name: str,
        category: str,
        cat: Dict[str, Any],
        warnings: List[WarningMessage],
    ) -> Optional[ResolvedFeature]:
        if not self._addgene.has_token():
            return None

        for kwargs in self._addgene_search_kwargs(clean_name, category):
            print(f"[FEATURE ADDGENE] '{clean_name}' search kwargs={kwargs}", flush=True)
            plasmids = self._addgene.search_plasmids(**kwargs)
            if not plasmids:
                continue

            for plasmid in plasmids[:3]:
                plasmid_id = plasmid.get("id")
                if not plasmid_id:
                    continue

                print(f"[FEATURE ADDGENE] fetching sequences for plasmid id={plasmid_id}", flush=True)
                plasmid_with_seq = self._addgene.get_plasmid_with_sequences(plasmid_id)
                if not plasmid_with_seq:
                    continue

                gb_text = self._addgene.download_first_full_genbank(plasmid_with_seq)
                if not gb_text:
                    continue

                seq = self._extract_seq(gb_text, clean_name, cat)
                if not seq:
                    continue

                label = f"Addgene:{plasmid_id}"
                if not self._validate_length(seq, clean_name, label, cat, warnings):
                    continue

                print(f"[FEATURE ADDGENE] ✓ '{clean_name}' plasmid={plasmid_id} ({len(seq)} bp)", flush=True)
                return ResolvedFeature(
                    name=clean_name,
                    type=cat["output_type"],
                    sequence=seq,
                    source=f"Addgene:{plasmid_id}",
                    length=len(seq),
                    validated=True,
                    position_hint=None,
                )

        return None

    # ------------------------------------------------------------------
    # Unified fetch: NCBI → Addgene → warn
    # ------------------------------------------------------------------

    def _fetch_feature(
        self,
        name: str,
        category: str,
        warnings: List[WarningMessage],
    ) -> Optional[ResolvedFeature]:
        clean_name = self._clean(name)
        cat = _CATEGORY[category]

        # 0. Known fixed accessions — bypass search entirely, raise on bad result
        feat = self._fetch_by_known_accession(clean_name, cat)
        if feat is not None:
            return feat

        # 1. NCBI dynamic search (primary)
        feat = self._try_ncbi(clean_name, category, cat, warnings)
        if feat:
            return feat

        # 2. Addgene (fallback)
        print(f"[FEATURE] '{clean_name}' not found on NCBI → trying Addgene", flush=True)
        feat = self._try_addgene(clean_name, category, cat, warnings)
        if feat:
            return feat

        warnings.append(WarningMessage(
            code="FEATURE_NOT_RESOLVED",
            message=f"Could not resolve '{name}' from NCBI or Addgene.",
        ))
        return None

    def _fetch_tag(
        self,
        name: str,
        position_hint: Optional[str],
        warnings: List[WarningMessage],
    ) -> Optional[ResolvedFeature]:
        feat = self._fetch_feature(name, "tag", warnings)
        if feat:
            return ResolvedFeature(
                name=feat.name,
                type=feat.type,
                sequence=feat.sequence,
                source=feat.source,
                length=feat.length,
                validated=feat.validated,
                position_hint=position_hint,
            )
        return None

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self, inp: FeatureInput) -> FeatureOutput:
        self._logger.info("[FEATURE] INPUT %s", inp.model_dump())
        warnings: List[WarningMessage] = []
        out_features: List[ResolvedFeature] = []

        host = (inp.expression_host or "").lower()
        is_mammalian = any(x in host for x in _MAMMALIAN_HOST_TOKENS)

        # 1. Kozak — fixed consensus, no database record needed
        if is_mammalian:
            out_features.append(ResolvedFeature(
                name="Kozak",
                type="kozak",
                sequence=KOZAK_SEQ,
                source="constant",
                length=len(KOZAK_SEQ),
                validated=True,
                position_hint=None,
            ))

        # 2. Promoter
        if inp.promoter:
            feat = self._fetch_feature(inp.promoter, "promoter", warnings)
            if feat:
                out_features.append(feat)

        # 3. Terminator / polyA (terminator field takes precedence)
        term_name = inp.terminator or inp.polyA
        category  = "terminator" if inp.terminator else "polya"
        if term_name:
            feat = self._fetch_feature(term_name, category, warnings)
            if feat:
                out_features.append(feat)
        elif is_mammalian:
            feat = self._fetch_feature("BGH polyA", "polya", warnings)
            if feat:
                out_features.append(feat)
                warnings.append(WarningMessage(
                    code="DEFAULT_TERMINATOR_ADDED",
                    message="No terminator/polyA provided. Added default BGH polyA for mammalian host.",
                ))

        # 4. Selection marker
        if inp.selection_marker:
            feat = self._fetch_feature(inp.selection_marker, "selection_marker", warnings)
            if feat:
                out_features.append(feat)

        # 5. Origin of replication
        if inp.origin_of_replication:
            feat = self._fetch_feature(inp.origin_of_replication, "ori", warnings)
            if feat:
                out_features.append(feat)

        # 6. Tags
        _pos_map: Dict[str, str] = {"N-terminal": "N", "C-terminal": "C"}
        tag_requests: List[Tuple[str, Optional[str]]] = []
        if inp.n_terminal_tag:
            tag_requests.append((inp.n_terminal_tag, "N"))
        if inp.c_terminal_tag:
            tag_requests.append((inp.c_terminal_tag, "C"))
        for ts in inp.extra_tags:
            tag_requests.append((ts.name, _pos_map.get(ts.position or "", None)))

        for tag_name, position in tag_requests:
            feat = self._fetch_tag(tag_name, position, warnings)
            if feat:
                out_features.append(feat)

        out = FeatureOutput(features=out_features, warnings=warnings)
        self._logger.info(
            "[FEATURE] OUTPUT features=%d warnings=%d",
            len(out.features), len(out.warnings),
        )
        return out
