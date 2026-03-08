from __future__ import annotations

import re
import urllib.parse
import urllib.request
from typing import Dict, Any, List, Optional, Tuple

from logging_utils import get_conversation_logger

from .models import BackboneInput, BackboneOutput, WarningMessage


class BackboneAgent:
    """
    NCBI-only public API backbone resolver.

    Workflow:
    1. Search NCBI nuccore for candidate vector/plasmid records.
    2. Fetch annotated GenBank flatfiles (gbwithparts).
    3. Parse FEATURES and sequence.
    4. Classify records from annotations:
       - clean backbone-like
       - loaded vector
    5. If loaded:
       - either reject and keep searching
       - or clean by removing the annotated expression-cassette/cargo span
    6. Return a real sequence only.

    Important:
    - No local library.
    - No synthetic fallback.
    - No title-only classification.
    """

    NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    NCBI_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    # Backbone-maintenance genes / elements that should not count as cargo
    SAFE_BACKBONE_MARKERS = {
        "bla", "amp", "ampr", "ampicillin", "betalactamase",
        "kan", "kana", "kanr", "kanamycin", "neo", "neor", "aph",
        "cat", "chloramphenicol", "tet", "tetracycline",
        "hygro", "hygromycin", "bsd", "blast", "zeo", "zeocin",
        "ori", "f1", "cole1", "puc", "pbr322", "rop",
        "mcs", "multiple cloning site", "polylinker",
        "lacz", "laczalpha", "rrn"
    }

    # Payload / cargo / engineered-expression markers
    PAYLOAD_MARKERS = {
        "gfp", "egfp", "sfgfp", "rfp", "cfp", "yfp", "bfp",
        "mcherry", "mscarlet", "mcardinal", "tdtomato", "venus",
        "reporter", "fluor", "fluorescent", "luciferase", "nanoluc",
        "fusion", "tag", "flag", "ha", "myc", "his", "strep", "v5",
        "cassette", "insert", "payload", "transgene",
        "nes", "nls", "devd", "tevp", "p2a", "t2a", "e2a", "ires",
        "cre", "cas9", "dcas9", "sensor", "biosensor",
        "tp53"
    }

    PROMOTER_MARKERS = {
        "cmv": "CMV",
        "ef1a": "EF1a",
        "ef1alpha": "EF1a",
        "cag": "CAG",
        "pgk": "PGK",
        "sv40": "SV40",
        "t7": "T7",
        "lac": "lac",
    }

    TERMINATOR_MARKERS = {
        "polya", "polyadenylation", "poly a", "bgh poly", "bgh polya",
        "sv40 poly", "sv40 polya", "terminator"
    }

    def __init__(
        self,
        *,
        ncbi_api_key: Optional[str] = None,
        ncbi_email: Optional[str] = None,
        user_agent: str = "plasmid-pipeline/0.1",
        allow_clean_loaded_vectors: bool = True,
    ) -> None:
        self._logger = get_conversation_logger()
        self._ncbi_api_key = ncbi_api_key
        self._ncbi_email = ncbi_email
        self._user_agent = user_agent
        self._allow_clean_loaded_vectors = allow_clean_loaded_vectors

    # ------------------------------------------------------------------
    # Basic helpers
    # ------------------------------------------------------------------

    def _clean_sequence(self, seq: Optional[str]) -> str:
        if not seq:
            return ""
        seq = seq.upper().replace("\n", "").replace("\r", "").replace(" ", "")
        return "".join(base for base in seq if base in {"A", "T", "G", "C", "N"})

    def _is_plausible_backbone_sequence(self, seq: str) -> bool:
        return len(self._clean_sequence(seq)) >= 1000

    def _normalize_text(self, text: str) -> str:
        t = (text or "").lower()
        for ch in [",", ";", "(", ")", "[", "]", "{", "}", "/", "|", ":", "-", "_"]:
            t = t.replace(ch, " ")
        t = re.sub(r"\s+", " ", t).strip()
        return t

    def _http_get_json(self, url: str, params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        query = urllib.parse.urlencode(params, doseq=True)
        full_url = f"{url}?{query}"
        req = urllib.request.Request(
            full_url,
            headers={"User-Agent": self._user_agent, "Accept": "application/json"},
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as response:
                return __import__("json").loads(response.read().decode("utf-8"))
        except Exception as e:
            self._logger.warning("[BACKBONE] JSON GET failed url=%s error=%s", full_url, str(e))
            return None

    def _http_get_text(self, url: str, params: Dict[str, Any]) -> Optional[str]:
        query = urllib.parse.urlencode(params, doseq=True)
        full_url = f"{url}?{query}"
        req = urllib.request.Request(
            full_url,
            headers={"User-Agent": self._user_agent, "Accept": "text/plain"},
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as response:
                return response.read().decode("utf-8")
        except Exception as e:
            self._logger.warning("[BACKBONE] TEXT GET failed url=%s error=%s", full_url, str(e))
            return None

    def _ncbi_common_params(self) -> Dict[str, Any]:
        params: Dict[str, Any] = {}
        if self._ncbi_api_key:
            params["api_key"] = self._ncbi_api_key
        if self._ncbi_email:
            params["email"] = self._ncbi_email
        return params

    # ------------------------------------------------------------------
    # Query construction
    # ------------------------------------------------------------------

    def _build_search_terms(self, inp: BackboneInput) -> List[str]:
        host = inp.expression_host.lower()
        promoter = (inp.promoter or "").upper()
        vector_type = (inp.vector_type or "").lower()

        terms: List[str] = []

        if ("hek" in host or "293" in host or "mamm" in host or "cho" in host):
            if promoter == "CMV":
                terms.extend([
                    'pcDNA[Title] AND vector[Title]',
                    '"mammalian expression vector"',
                    '"CMV expression vector"',
                    'plasmid[Title] AND mammalian[Title]',
                ])
            else:
                terms.extend([
                    '"mammalian expression vector"',
                    'pcDNA[Title] AND vector[Title]',
                    'plasmid[Title] AND mammalian[Title]',
                ])
        elif ("coli" in host or "bacteria" in host or "e. coli" in host):
            if promoter == "T7":
                terms.extend([
                    'pET[Title] AND vector[Title]',
                    '"T7 expression vector"',
                    '"bacterial expression vector"',
                ])
            else:
                terms.extend([
                    '"bacterial expression vector"',
                    'plasmid[Title] AND vector[Title]',
                ])
        else:
            terms.extend([
                '"expression vector"',
                'plasmid[Title] AND vector[Title]',
            ])

        if vector_type:
            terms.append(f'"{vector_type}"')

        terms.append('(plasmid[Title] OR vector[Title]) AND "complete sequence"')

        seen = set()
        out: List[str] = []
        for t in terms:
            if t not in seen:
                out.append(t)
                seen.add(t)
        return out

    # ------------------------------------------------------------------
    # NCBI retrieval
    # ------------------------------------------------------------------

    def _ncbi_esearch(self, term: str, retmax: int = 12) -> List[str]:
        data = self._http_get_json(
            self.NCBI_ESEARCH_URL,
            {
                "db": "nuccore",
                "term": term,
                "retmode": "json",
                "retmax": retmax,
                **self._ncbi_common_params(),
            },
        )
        if not data:
            return []
        return ((data.get("esearchresult") or {}).get("idlist")) or []

    def _ncbi_esummary(self, ids: List[str]) -> Dict[str, Any]:
        if not ids:
            return {}
        return self._http_get_json(
            self.NCBI_ESUMMARY_URL,
            {
                "db": "nuccore",
                "id": ",".join(ids),
                "retmode": "json",
                **self._ncbi_common_params(),
            },
        ) or {}

    def _extract_summary_docs(self, esummary: Dict[str, Any]) -> List[Tuple[str, Dict[str, Any]]]:
        result = esummary.get("result") or {}
        uids = result.get("uids") or []
        docs: List[Tuple[str, Dict[str, Any]]] = []
        for uid in uids:
            doc = result.get(uid)
            if isinstance(doc, dict):
                docs.append((uid, doc))
        return docs

    def _ncbi_fetch_genbank(self, nuccore_id: str) -> Optional[str]:
        return self._http_get_text(
            self.NCBI_EFETCH_URL,
            {
                "db": "nuccore",
                "id": nuccore_id,
                "rettype": "gbwithparts",
                "retmode": "text",
                **self._ncbi_common_params(),
            },
        )

    # ------------------------------------------------------------------
    # GenBank parsing
    # ------------------------------------------------------------------

    def _parse_genbank_sequence(self, gb_text: str) -> str:
        lines = gb_text.splitlines()
        in_origin = False
        seq_parts: List[str] = []

        for line in lines:
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin:
                if line.startswith("//"):
                    break
                seq_parts.append("".join(ch for ch in line if ch.isalpha()))

        return self._clean_sequence("".join(seq_parts))

    def _parse_genbank_features(self, gb_text: str) -> List[Dict[str, Any]]:
        lines = gb_text.splitlines()
        in_features = False
        features: List[Dict[str, Any]] = []
        current: Optional[Dict[str, Any]] = None
        current_qualifier: Optional[str] = None

        feature_key_re = re.compile(r"^\s{5}(\S+)\s+(.+)$")
        qualifier_re = re.compile(r'^\s+/([^=]+)(?:=(.*))?$')

        for line in lines:
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if in_features and line.startswith("ORIGIN"):
                break
            if not in_features:
                continue

            m = feature_key_re.match(line)
            if m:
                if current is not None:
                    features.append(current)
                current = {
                    "type": m.group(1).strip(),
                    "location": m.group(2).strip(),
                    "qualifiers": {},
                }
                current_qualifier = None
                continue

            if current is None:
                continue

            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith("/"):
                qm = qualifier_re.match(stripped)
                if not qm:
                    continue

                qkey = qm.group(1).strip()
                qval = qm.group(2)

                if qval is None:
                    value = ""
                else:
                    value = qval.strip()
                    if value.startswith('"') and value.endswith('"') and len(value) >= 2:
                        value = value[1:-1]

                current["qualifiers"].setdefault(qkey, []).append(value)
                current_qualifier = qkey
                continue

            if current_qualifier is not None:
                cont = stripped
                if cont.endswith('"'):
                    cont = cont[:-1]
                if cont.startswith('"'):
                    cont = cont[1:]

                existing = current["qualifiers"][current_qualifier]
                if existing:
                    existing[-1] = (existing[-1] + " " + cont).strip()
                else:
                    existing.append(cont)

        if current is not None:
            features.append(current)

        return features

    def _feature_text(self, feature: Dict[str, Any]) -> str:
        quals = feature.get("qualifiers", {})
        parts: List[str] = [feature.get("type", "")]
        for qvals in quals.values():
            for v in qvals:
                parts.append(v)
        return self._normalize_text(" ".join(parts))

    def _feature_hits(self, text: str, vocab: set[str]) -> List[str]:
        hits: List[str] = []
        for marker in vocab:
            if marker in text:
                hits.append(marker)
        return sorted(set(hits))

    def _location_bounds(self, location: str) -> Optional[Tuple[int, int]]:
        """
        Parse a GenBank location string into (start, end), 1-based inclusive.
        Works for simple/complement/join-style numeric locations by taking min/max coordinates.
        """
        nums = [int(x) for x in re.findall(r"\d+", location or "")]
        if len(nums) < 2:
            return None
        return min(nums), max(nums)

    def _extract_promoters_from_features(self, features: List[Dict[str, Any]]) -> List[str]:
        promoters: List[str] = []

        for feat in features:
            text = self._feature_text(feat)
            if feat.get("type") == "promoter" or "promoter" in text:
                for k, label in self.PROMOTER_MARKERS.items():
                    if k in text and label not in promoters:
                        promoters.append(label)

        for feat in features:
            text = self._feature_text(feat)
            for k, label in self.PROMOTER_MARKERS.items():
                if k in text and label not in promoters:
                    promoters.append(label)

        return promoters

    def _extract_terminator_like_features(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        out: List[Dict[str, Any]] = []
        for feat in features:
            text = self._feature_text(feat)
            if feat.get("type") in {"terminator", "polyA_signal", "regulatory"}:
                out.append(feat)
                continue
            if any(marker in text for marker in self.TERMINATOR_MARKERS):
                out.append(feat)
        return out

    # ------------------------------------------------------------------
    # Annotation-based classification
    # ------------------------------------------------------------------

    def _classify_record(self, gb_text: str) -> Dict[str, Any]:
        features = self._parse_genbank_features(gb_text)
        promoters = self._extract_promoters_from_features(features)

        cds_total = 0
        non_backbone_cds = 0
        payload_markers: List[str] = []

        for feat in features:
            ftype = feat.get("type", "")
            text = self._feature_text(feat)

            if ftype == "CDS":
                cds_total += 1
                safe_hits = self._feature_hits(text, self.SAFE_BACKBONE_MARKERS)
                payload_hits = self._feature_hits(text, self.PAYLOAD_MARKERS)
                if payload_hits:
                    payload_markers.extend(payload_hits)
                if not safe_hits:
                    non_backbone_cds += 1

            elif ftype in {"gene", "misc_feature", "regulatory", "promoter"}:
                payload_hits = self._feature_hits(text, self.PAYLOAD_MARKERS)
                if payload_hits:
                    payload_markers.extend(payload_hits)

        payload_markers = sorted(set(payload_markers))

        is_loaded = False
        if payload_markers:
            is_loaded = True
        elif non_backbone_cds > 1:
            is_loaded = True
        elif non_backbone_cds == 1 and promoters:
            is_loaded = True

        return {
            "features": features,
            "promoters": promoters,
            "payload_markers": payload_markers,
            "cds_total": cds_total,
            "non_backbone_cds": non_backbone_cds,
            "is_loaded": is_loaded,
        }

    # ------------------------------------------------------------------
    # Cleaning loaded vectors from annotation
    # ------------------------------------------------------------------

    def _find_cargo_span(self, parsed: Dict[str, Any]) -> Optional[Tuple[int, int]]:
        """
        Try to identify the expression-cassette / cargo span from annotations only.

        Strategy:
        - collect suspicious features:
          promoter + non-backbone CDS + payload-marked misc/gene/regulatory + polyA/terminator
        - if they form a coherent annotated block, remove the full span
        """
        features = parsed["features"]
        suspicious_bounds: List[Tuple[int, int]] = []

        terminator_like = self._extract_terminator_like_features(features)

        for feat in features:
            text = self._feature_text(feat)
            bounds = self._location_bounds(feat.get("location", ""))
            if bounds is None:
                continue

            ftype = feat.get("type", "")
            safe_hits = self._feature_hits(text, self.SAFE_BACKBONE_MARKERS)
            payload_hits = self._feature_hits(text, self.PAYLOAD_MARKERS)

            add = False

            if ftype == "promoter":
                add = True
            elif "promoter" in text and ftype in {"misc_feature", "regulatory"}:
                add = True
            elif ftype == "CDS" and not safe_hits:
                add = True
            elif payload_hits:
                add = True

            if add:
                suspicious_bounds.append(bounds)

        for feat in terminator_like:
            bounds = self._location_bounds(feat.get("location", ""))
            if bounds is not None:
                suspicious_bounds.append(bounds)

        if not suspicious_bounds:
            return None

        start = min(b[0] for b in suspicious_bounds)
        end = max(b[1] for b in suspicious_bounds)

        if end <= start:
            return None

        return start, end

    def _clean_loaded_backbone(
        self,
        seq: str,
        parsed: Dict[str, Any],
    ) -> Tuple[Optional[str], List[WarningMessage], Optional[int]]:
        warnings: List[WarningMessage] = []

        cargo_span = self._find_cargo_span(parsed)
        if cargo_span is None:
            warnings.append(
                WarningMessage(
                    code="LOADED_BACKBONE_CANNOT_BE_CLEANED",
                    message="The record looks loaded, but no coherent annotated cargo span could be identified for removal.",
                )
            )
            return None, warnings, None

        start_1, end_1 = cargo_span

        # Convert 1-based inclusive -> python slicing
        start_0 = max(0, start_1 - 1)
        end_0 = min(len(seq), end_1)

        left = seq[:start_0]
        right = seq[end_0:]
        cleaned = left + right

        if not self._is_plausible_backbone_sequence(cleaned):
            warnings.append(
                WarningMessage(
                    code="CLEANED_BACKBONE_TOO_SHORT",
                    message=(
                        f"Annotated cargo removal produced a sequence of length {len(cleaned)}, "
                        "which is too short to trust as a real backbone."
                    ),
                )
            )
            return None, warnings, None

        warnings.append(
            WarningMessage(
                code="LOADED_BACKBONE_CLEANED",
                message=(
                    f"Loaded vector was converted into a putative bare backbone by removing "
                    f"annotated cargo span {start_1}..{end_1}."
                ),
            )
        )

        insertion_site_index = len(left)
        return cleaned, warnings, insertion_site_index

    # ------------------------------------------------------------------
    # Candidate scoring and selection
    # ------------------------------------------------------------------

    def _score_candidate_from_annotations(
        self,
        inp: BackboneInput,
        summary_doc: Dict[str, Any],
        parsed: Dict[str, Any],
        seq_len: int,
    ) -> int:
        score = 0

        title = self._normalize_text(str(summary_doc.get("title") or ""))
        host = inp.expression_host.lower()
        promoter = (inp.promoter or "").upper()

        if seq_len >= 1000:
            score += 20
        if seq_len >= 3000:
            score += 10

        if "vector" in title:
            score += 10
        if "plasmid" in title:
            score += 10
        if "cloning" in title:
            score += 8
        if "complete sequence" in title:
            score += 2

        promoters = parsed["promoters"]
        if promoter and promoter in promoters:
            score += 15

        if "hek" in host or "293" in host or "mamm" in host or "cho" in host:
            if "pcdna" in title:
                score += 10
            if "mammalian" in title:
                score += 8
        elif "coli" in host or "bacteria" in host:
            if "pet" in title:
                score += 10
            if "bacterial" in title:
                score += 8

        # Strong penalty to loaded vectors as raw candidates
        if parsed["is_loaded"]:
            score -= 100
        score -= 20 * parsed["non_backbone_cds"]
        score -= 5 * len(parsed["payload_markers"])

        return score

    def _fetch_remote_backbone(self, inp: BackboneInput) -> Tuple[Optional[Dict[str, Any]], List[WarningMessage]]:
        warnings: List[WarningMessage] = []
        terms = self._build_search_terms(inp)

        best_clean_candidate: Optional[Dict[str, Any]] = None
        best_clean_score = -10**9

        best_cleaned_loaded_candidate: Optional[Dict[str, Any]] = None
        best_cleaned_loaded_score = -10**9

        for term in terms:
            ids = self._ncbi_esearch(term, retmax=12)
            if not ids:
                continue

            esummary = self._ncbi_esummary(ids)
            docs = self._extract_summary_docs(esummary)

            for uid, doc in docs:
                gb_text = self._ncbi_fetch_genbank(uid)
                if not gb_text:
                    continue

                seq = self._parse_genbank_sequence(gb_text)
                if not self._is_plausible_backbone_sequence(seq):
                    continue

                parsed = self._classify_record(gb_text)
                score = self._score_candidate_from_annotations(
                    inp=inp,
                    summary_doc=doc,
                    parsed=parsed,
                    seq_len=len(seq),
                )

                base_notes = [
                    f"Resolved from NCBI nuccore record {uid}.",
                    f"Search term used: {term}",
                    f"Annotation-derived promoters: {parsed['promoters']}",
                    f"Annotation-derived payload markers: {parsed['payload_markers']}",
                    f"Annotation-derived non-backbone CDS count: {parsed['non_backbone_cds']}",
                ]

                # Case 1: already clean enough
                if not parsed["is_loaded"]:
                    if score > best_clean_score:
                        best_clean_candidate = {
                            "name": str(doc.get("title") or f"NCBI_nuccore_{uid}"),
                            "sequence": seq,
                            "source": "ncbi_nuccore",
                            "accession": str(doc.get("caption") or uid),
                            "notes": base_notes,
                            "is_loaded_vector": False,
                            "backbone_promoters": parsed["promoters"],
                            "backbone_payload_markers": parsed["payload_markers"],
                            "suggested_strategy": "full_cassette_ok",
                            "insertion_site_index": len(seq) // 2,
                        }
                        best_clean_score = score
                    continue

                # Case 2: loaded, try to clean if allowed
                if self._allow_clean_loaded_vectors:
                    cleaned_seq, clean_warnings, insertion_site_index = self._clean_loaded_backbone(seq, parsed)
                    if cleaned_seq and self._is_plausible_backbone_sequence(cleaned_seq):
                        cleaned_score = score + 80  # recovering from a loaded record, but still below a naturally clean one unless annotations are good
                        candidate_notes = list(base_notes)
                        candidate_notes.append("Loaded record was cleaned using annotated GenBank feature span removal.")

                        candidate = {
                            "name": str(doc.get("title") or f"NCBI_nuccore_{uid}") + " [cleaned]",
                            "sequence": cleaned_seq,
                            "source": "ncbi_nuccore_cleaned",
                            "accession": str(doc.get("caption") or uid),
                            "notes": candidate_notes,
                            "is_loaded_vector": False,
                            "backbone_promoters": [],
                            "backbone_payload_markers": [],
                            "suggested_strategy": "full_cassette_ok_after_cleaning",
                            "insertion_site_index": insertion_site_index if insertion_site_index is not None else len(cleaned_seq) // 2,
                            "extra_warnings": clean_warnings,
                        }
                        if cleaned_score > best_cleaned_loaded_score:
                            best_cleaned_loaded_candidate = candidate
                            best_cleaned_loaded_score = cleaned_score

        chosen: Optional[Dict[str, Any]] = None

        if best_clean_candidate is not None:
            chosen = best_clean_candidate
        elif best_cleaned_loaded_candidate is not None:
            chosen = best_cleaned_loaded_candidate

        if chosen is None:
            warnings.append(
                WarningMessage(
                    code="NCBI_BACKBONE_NOT_FOUND",
                    message="No suitable clean backbone could be resolved from NCBI, and no loaded candidate could be safely cleaned.",
                )
            )
            return None, warnings

        warnings.append(
            WarningMessage(
                code="NCBI_BACKBONE_SUCCESS",
                message=f"Backbone resolved from NCBI: {chosen['name']}",
            )
        )

        for w in chosen.get("extra_warnings", []):
            warnings.append(w)

        return chosen, warnings

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self, inp: BackboneInput) -> BackboneOutput:
        self._logger.info("[BACKBONE] INPUT %s", inp.model_dump())

        chosen, warnings = self._fetch_remote_backbone(inp)

        if chosen is None:
            raise ValueError(
                "No validated clean backbone sequence could be resolved from public APIs."
            )

        seq = self._clean_sequence(chosen.get("sequence"))
        if not self._is_plausible_backbone_sequence(seq):
            raise ValueError(
                f"Resolved backbone '{chosen.get('name', 'unknown')}' is too short ({len(seq)} bp) to be treated as a real plasmid backbone."
            )

        notes = list(chosen.get("notes", []))
        if inp.promoter:
            notes.append(f"Requested promoter={inp.promoter}.")
        if inp.expression_host:
            notes.append(f"Requested expression_host={inp.expression_host}.")
        if inp.vector_type:
            notes.append(f"Requested vector_type={inp.vector_type}.")

        out = BackboneOutput(
            backbone_name=str(chosen.get("name") or "resolved_backbone"),
            backbone_sequence=seq,
            source=str(chosen.get("source") or "unknown"),
            notes=notes,
            warnings=warnings,
            is_loaded_vector=bool(chosen.get("is_loaded_vector", False)),
            backbone_promoters=list(chosen.get("backbone_promoters", [])),
            backbone_payload_markers=list(chosen.get("backbone_payload_markers", [])),
            suggested_strategy=chosen.get("suggested_strategy"),
            insertion_site_index=chosen.get("insertion_site_index"),
        )

        self._logger.info(
            "[BACKBONE] OUTPUT name=%s source=%s length=%s loaded=%s warnings=%s",
            out.backbone_name,
            out.source,
            len(out.backbone_sequence),
            out.is_loaded_vector,
            [w.model_dump() for w in out.warnings],
        )
        return out