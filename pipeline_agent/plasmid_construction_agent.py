from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import httpx

from logging_utils import get_conversation_logger
from .models import (
    IntentOutput,
    PlasmidConstructionInput,
    PlasmidConstructionOutput,
    ResolvedElement,
    WEB_FALLBACK_WARNING,
)


@dataclass(frozen=True)
class _Slot:
    name: str
    element_type: str
    requested_name: Optional[str]


class PlasmidConstructionAgent:
    """
    Builds a plasmid in template order using:
    1) backbone retrieval & feature detection
    2) local dataset lookup
    3) NCBI web fallback for missing elements.
    """

    _DNA_COLUMNS = (
        "DNA Sequence (5' → 3')",
        "DNA Sequence",
        "sequence",
    )
    _AA_COLUMNS = (
        "Amino Acid Sequence",
        "Tag 1 (aa)",
        "Tag 2 (aa)",
        "Linker (aa)",
    )
    _CDS_ELEMENT_TYPES = ("CDS", "gene", "cargo", "reporter", "insert")

    # NCBI sequence-length filters per element type to avoid fetching whole plasmids
    # when searching for sub-elements.  Format: "min:max" passed to [Sequence Length].
    _ELEMENT_SIZE_RANGES: Dict[str, str] = {
        "rep_origin":    "400:5000",
        "polyA_signal":  "150:700",
        "terminator":    "20:600",
        "promoter":      "20:2000",
        "misc_feature":  "10:2000",
        "RBS":           "5:100",
    }
    # Minimum backbone size — avoids picking up short partial MCS records.
    _BACKBONE_MIN_SIZE = 2000

    # Columns that may carry a direct NCBI accession or Addgene plasmid ID.
    _ACCESSION_COLUMNS = ("NCBI Accession", "GenBank Accession", "Accession")
    _ADDGENE_COLUMNS   = ("Addgene ID", "Addgene Plasmid ID")

    _NAME_COLUMNS = (
        "Element Name",
        "Backbone Name",
        "Tag Name",
        "Linker Name",
        "Terminator Name",
        "Fusion Name",
        "NLS Type",
        "Cargo Type",
        "Cargo Name",
        "Reporter Name",
    )
    _TYPE_COLUMNS = ("Type", "Origin Type")
    _CODON_TABLE = {
        "A": "GCT",
        "C": "TGT",
        "D": "GAT",
        "E": "GAA",
        "F": "TTT",
        "G": "GGT",
        "H": "CAT",
        "I": "ATT",
        "K": "AAA",
        "L": "CTG",
        "M": "ATG",
        "N": "AAT",
        "P": "CCT",
        "Q": "CAA",
        "R": "CGT",
        "S": "TCT",
        "T": "ACT",
        "V": "GTG",
        "W": "TGG",
        "Y": "TAT",
        "*": "TAA",
    }

    def __init__(
        self,
        *,
        dataset_csv: Optional[str] = None,
        timeout_s: float = 20.0,
        strict: bool = False,
    ) -> None:
        self._logger = get_conversation_logger()
        self._timeout_s = timeout_s
        root = Path(__file__).resolve().parents[1]
        self._dataset_dir = str(root / "data")
        self._dataset_csv = dataset_csv
        self._dataset = self._load_dataset(dataset_csv)
        self._addgene_api_key = os.getenv("ADDGENE_API_KEY")
        # Strict mode: never emit placeholder sequences.  Controlled by constructor
        # param or PLASMID_STRICT_MODE env var ("1" / "true" / "yes").
        _env = os.getenv("PLASMID_STRICT_MODE", "").strip().lower()
        self._strict = strict or (_env in ("1", "true", "yes"))

    @staticmethod
    def _sanitize_dna(seq: str) -> str:
        cleaned = re.sub(r"[^ACGTNacgtn]", "", seq).upper()
        return cleaned

    @staticmethod
    def _placeholder_seq(name: str, element_type: str) -> str:
        seed = (name + element_type).upper()
        return ("NNNN" + re.sub(r"[^ACGT]", "", seed) + "NNNN")[:120] or ("N" * 60)

    @staticmethod
    def _sanitize_text(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", " ", (value or "").lower()).strip()

    def _load_dataset(self, single_csv: Optional[str]) -> List[Dict[str, str]]:
        paths: List[Path] = []
        if single_csv:
            p = Path(single_csv)
            if p.exists():
                paths = [p]
        else:
            data_dir = Path(self._dataset_dir)
            if data_dir.exists():
                paths = sorted(data_dir.glob("plasmid_*.csv"))

        if not paths:
            self._logger.info("[CONSTRUCTION] no dataset csv files found.")
            return []

        rows: List[Dict[str, str]] = []
        for p in paths:
            with open(p, "r", encoding="utf-8", newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    item = dict(row)
                    item["_dataset_file"] = p.name
                    rows.append(item)
        self._logger.info("[CONSTRUCTION] loaded %d rows from %d csv files", len(rows), len(paths))
        return rows

    @staticmethod
    def _first_nonempty(row: Dict[str, str], keys: Tuple[str, ...]) -> str:
        for key in keys:
            val = (row.get(key) or "").strip()
            if val:
                return val
        return ""

    def _get_row_name(self, row: Dict[str, str]) -> str:
        return self._first_nonempty(row, self._NAME_COLUMNS) or self._first_nonempty(
            row, tuple(k for k in row.keys() if not k.startswith("_"))
        )

    def _extract_dna(self, row: Dict[str, str]) -> str:
        raw = self._first_nonempty(row, self._DNA_COLUMNS)
        return self._sanitize_dna(raw)

    @staticmethod
    def _normalize_aa(value: str) -> str:
        if not value:
            return ""
        seq = value.strip().upper()
        seq = re.sub(r"\[.*?\]", "", seq)
        seq = seq.replace("...", "")
        # Expand common repeat notation like (GGGGS)n to 3 repeats.
        m = re.fullmatch(r"\(([-A-Z*]+)\)N", seq)
        if m:
            seq = m.group(1) * 3
        return re.sub(r"[^A-Z*]", "", seq)

    def _translate_aa_to_dna(self, aa_seq: str) -> str:
        aa = self._normalize_aa(aa_seq)
        if not aa or any(ch not in self._CODON_TABLE for ch in aa):
            return ""
        return "".join(self._CODON_TABLE[ch] for ch in aa)

    def _extract_or_build_dna_from_row(self, row: Dict[str, str]) -> str:
        dna = self._extract_dna(row)
        if dna:
            return dna
        aa_raw = self._first_nonempty(row, self._AA_COLUMNS)
        return self._translate_aa_to_dna(aa_raw)

    def _get_row_accession(self, row: Dict[str, str]) -> Tuple[str, str]:
        """Return (accession_value, source_type) from a row, or ('', '')."""
        for col in self._ACCESSION_COLUMNS:
            val = (row.get(col) or "").strip()
            if val:
                return val, "ncbi"
        for col in self._ADDGENE_COLUMNS:
            val = (row.get(col) or "").strip()
            if val:
                return val, "addgene"
        return "", ""

    def _match_dataset(
        self, requested_name: str, element_type: str
    ) -> Tuple[Optional[ResolvedElement], str, str]:
        """
        Returns (resolved_element, accession_hint, accession_type).

        resolved_element  – a ResolvedElement with sequence if dataset row has DNA/AA.
        accession_hint    – NCBI accession or Addgene ID from the best matching row
                            that had no sequenceable data.  Empty string when unused.
        accession_type    – 'ncbi' | 'addgene' | ''.
        """
        q = self._sanitize_text(requested_name)
        t = element_type.strip().lower()
        best: Optional[Tuple[int, ResolvedElement]] = None
        best_acc: str = ""
        best_acc_type: str = ""
        best_acc_score: int = -1

        for row in self._dataset:
            feat = self._get_row_name(row)
            typ = self._first_nonempty(row, self._TYPE_COLUMNS)
            if not feat:
                continue
            feat_l = self._sanitize_text(feat)
            typ_l = typ.lower()
            if not (q == feat_l or q in feat_l or feat_l in q):
                continue

            score = 3 if q == feat_l else 2
            if t and typ_l and (t in typ_l or typ_l in t):
                score += 1

            seq = self._extract_or_build_dna_from_row(row)
            if seq:
                candidate = ResolvedElement(
                    slot=element_type,
                    requested_name=requested_name,
                    element_type=typ or element_type,
                    sequence=seq,
                    source=f"dataset:{row.get('_dataset_file', 'unknown.csv')}",
                    from_web=False,
                )
                if best is None or score > best[0]:
                    best = (score, candidate)
            else:
                # Row has metadata / accession but no inline sequence.
                acc, acc_type = self._get_row_accession(row)
                if acc and score > best_acc_score:
                    best_acc = acc
                    best_acc_type = acc_type
                    best_acc_score = score

        resolved = best[1] if best else None
        # Prefer the accession from the sequence-bearing row itself when available.
        if resolved:
            row_acc, row_acc_type = "", ""
            # We don't store the row reference, so skip; accession hint is only needed
            # when resolved is None.
            return resolved, "", ""
        return None, best_acc, best_acc_type

    def _search_addgene_sequence(self, requested_name: str) -> Tuple[str, str]:
        if not self._addgene_api_key:
            raise ValueError("ADDGENE_API_KEY not set.")
        headers = {"Authorization": f"Token {self._addgene_api_key}"}
        with httpx.Client(
            timeout=self._timeout_s,
            follow_redirects=True,
            base_url="https://api.developers.addgene.org",
            headers=headers,
        ) as client:
            search = client.get("/catalog/plasmid/", params={"name": requested_name, "page_size": 1})
            search.raise_for_status()
            results = search.json().get("results") or []
            if not results:
                raise ValueError(f"No Addgene hit for '{requested_name}'.")
            plasmid_id = results[0].get("id")
            if not plasmid_id:
                raise ValueError("Addgene result has no plasmid id.")

            detail = client.get(f"/catalog/plasmid-with-sequences/{plasmid_id}/")
            detail.raise_for_status()
            payload = detail.json()
            seq_groups = payload.get("sequences", {})
            for key in (
                "public_user_full_sequences",
                "public_addgene_full_sequences",
                "public_user_partial_sequences",
                "public_addgene_partial_sequences",
            ):
                for item in (seq_groups.get(key) or []):
                    seq = self._sanitize_dna(item.get("sequence", ""))
                    if seq:
                        return seq, f"web:addgene:{plasmid_id}"
        raise ValueError(f"No Addgene sequence body for '{requested_name}'.")

    # ---------------------------------------------------------------------------
    # Addgene backbone query — sequence + annotated features
    # ---------------------------------------------------------------------------

    def _query_addgene_backbone(self, backbone_name: str) -> Optional[Dict]:
        """
        Fetch backbone sequence and feature annotations from Addgene.

        Tries two sources in order:
        1. Addgene developer API (ADDGENE_API_KEY required) — structured JSON
           with annotated features embedded in the sequence records.
        2. Addgene public vector database — no key required, feature data varies.

        Returns:
            {id, name, sequence, features: [{type, name, start, end}]}
        or None if both sources fail.
        """
        if self._addgene_api_key:
            result = self._query_addgene_developer_api(backbone_name)
            if result and result.get("sequence"):
                return result

        return self._query_addgene_vector_db(backbone_name)

    def _query_addgene_developer_api(self, backbone_name: str) -> Optional[Dict]:
        """Fetch backbone info via Addgene developer API with ADDGENE_API_KEY."""
        if not self._addgene_api_key:
            return None
        headers = {"Authorization": f"Token {self._addgene_api_key}"}
        try:
            with httpx.Client(
                timeout=self._timeout_s,
                follow_redirects=True,
                base_url="https://api.developers.addgene.org",
                headers=headers,
            ) as client:
                search_resp = client.get(
                    "/catalog/plasmid/",
                    params={"name": backbone_name, "page_size": 3},
                )
                search_resp.raise_for_status()
                results = search_resp.json().get("results") or []
                if not results:
                    self._logger.info(
                        "[CONSTRUCTION] Addgene API: no hit for backbone '%s'", backbone_name
                    )
                    return None
                plasmid_id = results[0].get("id")
                if not plasmid_id:
                    return None

                detail_resp = client.get(f"/catalog/plasmid-with-sequences/{plasmid_id}/")
                detail_resp.raise_for_status()
                payload = detail_resp.json()

                seq, features = self._unpack_addgene_payload(payload)
                self._logger.info(
                    "[CONSTRUCTION] Addgene API: '%s' → id=%s, seq=%d bp, features=%d",
                    backbone_name, plasmid_id, len(seq), len(features),
                )
                return {
                    "id": plasmid_id,
                    "name": payload.get("name", backbone_name),
                    "sequence": seq,
                    "features": features,
                }
        except Exception as exc:
            self._logger.info(
                "[CONSTRUCTION] Addgene developer API failed for '%s': %s", backbone_name, exc
            )
            return None

    def _query_addgene_vector_db(self, backbone_name: str) -> Optional[Dict]:
        """
        Query Addgene public vector database.
        Endpoint: https://www.addgene.org/vector-database/query/?q_vdb="{name}"
        Used as a fallback when no developer API key is available.
        """
        try:
            with httpx.Client(timeout=self._timeout_s, follow_redirects=True) as client:
                resp = client.get(
                    "https://www.addgene.org/vector-database/query/",
                    params={"q_vdb": f'"{backbone_name}"'},
                    headers={"Accept": "application/json, text/html;q=0.5"},
                )
                resp.raise_for_status()
                ct = resp.headers.get("content-type", "")
                if "json" not in ct:
                    self._logger.info(
                        "[CONSTRUCTION] Addgene vector DB returned non-JSON for '%s'", backbone_name
                    )
                    return None
                data = resp.json()
                items = (
                    data if isinstance(data, list)
                    else data.get("results") or data.get("vectors") or []
                )
                if not items:
                    return None
                item = items[0]
                seq = self._sanitize_dna(
                    item.get("sequence") or item.get("full_sequence") or ""
                )
                features = item.get("features") or item.get("annotations") or []
                self._logger.info(
                    "[CONSTRUCTION] Addgene vector DB: '%s' → seq=%d bp, features=%d",
                    backbone_name, len(seq), len(features),
                )
                return {
                    "id": item.get("id") or item.get("addgene_id"),
                    "name": item.get("name") or backbone_name,
                    "sequence": seq,
                    "features": features,
                }
        except Exception as exc:
            self._logger.info(
                "[CONSTRUCTION] Addgene vector DB failed for '%s': %s", backbone_name, exc
            )
            return None

    def _unpack_addgene_payload(self, payload: Dict) -> Tuple[str, List[Dict]]:
        """
        Extract the longest DNA sequence and feature list from an Addgene
        /catalog/plasmid-with-sequences/ response dict.

        Handles three sub-formats:
          • GenBank text in the "sequence" key (parsed with BioPython)
          • Inline annotations key ("annotations" / "features") in each record
          • Plain FASTA-like DNA string
        """
        best_seq = ""
        all_features: List[Dict] = []

        seq_groups = payload.get("sequences", {})
        for key in (
            "public_addgene_full_sequences",
            "public_user_full_sequences",
            "public_addgene_partial_sequences",
            "public_user_partial_sequences",
        ):
            for item in (seq_groups.get(key) or []):
                raw = (item.get("sequence") or "").strip()
                # GenBank text starts with LOCUS
                if raw.startswith("LOCUS"):
                    parsed_feats = self._parse_genbank_features(raw)
                    if parsed_feats and not all_features:
                        all_features = parsed_feats
                    dna = self._sanitize_dna(raw)
                else:
                    dna = self._sanitize_dna(raw)

                if len(dna) > len(best_seq):
                    best_seq = dna

                # Collect any inline annotation keys
                for feat in (item.get("annotations") or item.get("features") or []):
                    all_features.append(feat)

        return best_seq, all_features

    @staticmethod
    def _parse_genbank_features(gbk_text: str) -> List[Dict]:
        """
        Parse a GenBank flat-file string and return annotated features as dicts.
        Uses BioPython; returns [] if unavailable or parse fails.
        """
        if not gbk_text or "LOCUS" not in gbk_text[:200]:
            return []
        try:
            from io import StringIO
            from Bio import SeqIO
            record = SeqIO.read(StringIO(gbk_text), "genbank")
            out: List[Dict] = []
            for feat in record.features:
                if feat.type in ("source",):
                    continue
                label = (
                    (feat.qualifiers.get("label") or [""])[0]
                    or (feat.qualifiers.get("gene") or [""])[0]
                    or (feat.qualifiers.get("product") or [""])[0]
                    or (feat.qualifiers.get("note") or [""])[0]
                )
                out.append({
                    "type": feat.type,
                    "name": label,
                    "start": int(feat.location.start),
                    "end": int(feat.location.end),
                })
            return out
        except Exception:
            return []

    def _enrich_features_from_annotations(
        self,
        base_features: Dict,
        api_features: List[Dict],
    ) -> Dict:
        """
        Overlay structured Addgene/GenBank feature annotations on top of the
        sequence-scan baseline returned by _detect_features().

        Populates promoter_names / terminator_names / marker_names / ori_names
        and updates the boolean flags accordingly.
        """
        if not api_features:
            return base_features

        promoter_names: List[str] = []
        terminator_names: List[str] = []
        marker_names: List[str] = []
        ori_names: List[str] = []

        _MARKER_KW = (
            "ampicillin", "ampr", "bla",
            "kanamycin", "kanr",
            "neomycin", "hygromycin", "puromycin",
            "chloramphenicol", "cmr",
            "spectinomycin", "blasticidin", "zeocin",
            "tetracycline", "tetr", "g418",
        )
        _MAMMALIAN_MARKER_KW = (
            "neomycin", "hygromycin", "puromycin", "blasticidin", "zeocin", "g418",
        )

        for feat in api_features:
            ft = (feat.get("type") or "").lower().replace("-", " ").replace("_", " ")
            fn = (feat.get("name") or feat.get("label") or "").strip()
            fn_lower = fn.lower()

            if "promoter" in ft:
                if fn:
                    promoter_names.append(fn)
            elif "terminat" in ft or "polya" in ft or "polyadenyl" in ft:
                if fn:
                    terminator_names.append(fn)
            elif "origin" in ft or ft in ("rep origin", "ori"):
                if fn:
                    ori_names.append(fn)
            elif ft in ("cds", "gene", "misc feature"):
                if any(kw in fn_lower for kw in _MARKER_KW):
                    marker_names.append(fn)

        mammalian_marker = any(
            any(kw in m.lower() for kw in _MAMMALIAN_MARKER_KW)
            for m in marker_names
        )
        mammalian_promoter = base_features.get("has_mammalian_promoter") or any(
            any(kw in n.lower() for kw in ("cmv", "sv40", "cag", "ef1", "pgk", "ubc", "rous"))
            for n in promoter_names
        )

        enriched = dict(base_features)
        enriched.update({
            "has_promoter":          enriched.get("has_promoter") or bool(promoter_names),
            "has_mammalian_promoter": mammalian_promoter,
            "has_polyA":             enriched.get("has_polyA") or bool(terminator_names),
            "has_resistance_marker": enriched.get("has_resistance_marker") or bool(marker_names),
            "has_mammalian_marker":  enriched.get("has_mammalian_marker") or mammalian_marker,
            "promoter_names":   promoter_names,
            "terminator_names": terminator_names,
            "marker_names":     marker_names,
            "ori_names":        ori_names,
            "feature_source": "addgene_api",
        })
        return enriched

    @staticmethod
    def _parse_fasta_all(fasta_text: str) -> List[str]:
        """Parse a multi-record FASTA and return a list of raw sequences (one per record)."""
        seqs: List[str] = []
        current: List[str] = []
        for line in fasta_text.splitlines():
            line = line.strip()
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                    current = []
            elif line:
                current.append(line)
        if current:
            seqs.append("".join(current))
        return seqs

    def _search_web_cds(self, gene_name: str) -> Tuple[str, str]:
        """
        Fetch a clean CDS nucleotide sequence from NCBI for a given gene name.

        Strategy:
        1. Search nuccore for a RefSeq mRNA or gene record matching the name.
        2. Fetch with rettype=fasta_cds_nt to extract annotated CDS features only.
        3. Take the first (shortest valid) CDS — avoids returning whole-plasmid sequences.
        4. Fall back to plain fasta of the record if no annotated CDS found.
        """
        # Prefer RefSeq mRNA records (NM_*) — these are clean single-gene sequences.
        term = (
            f"{gene_name}[Gene Name] AND refseq[filter] AND 300:10000[Sequence Length]"
        )
        with httpx.Client(timeout=self._timeout_s, follow_redirects=True) as client:
            esearch = client.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params={"db": "nuccore", "retmode": "json", "retmax": 3, "term": term},
            )
            esearch.raise_for_status()
            ids = esearch.json().get("esearchresult", {}).get("idlist", [])

            # Fallback: broader search if RefSeq-specific gave nothing
            if not ids:
                term2 = f"{gene_name}[Title] AND 300:10000[Sequence Length]"
                esearch2 = client.get(
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                    params={"db": "nuccore", "retmode": "json", "retmax": 3, "term": term2},
                )
                esearch2.raise_for_status()
                ids = esearch2.json().get("esearchresult", {}).get("idlist", [])

            if not ids:
                raise ValueError(f"No NCBI hit for gene '{gene_name}'.")

            nuccore_id = ids[0]

            # Try fasta_cds_nt — returns individual FASTA per annotated CDS feature.
            efetch_cds = client.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params={
                    "db": "nuccore",
                    "id": nuccore_id,
                    "rettype": "fasta_cds_nt",
                    "retmode": "text",
                },
            )
            efetch_cds.raise_for_status()
            cds_records = self._parse_fasta_all(efetch_cds.text)
            for raw in cds_records:
                seq = self._sanitize_dna(raw)
                if len(seq) >= 300:  # skip tiny fragments
                    return seq, f"web:ncbi:{nuccore_id}:cds"

            # Fallback: plain FASTA of the record (e.g. short synthetic sequences)
            efetch_plain = client.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params={
                    "db": "nuccore",
                    "id": nuccore_id,
                    "rettype": "fasta",
                    "retmode": "text",
                },
            )
            efetch_plain.raise_for_status()
            plain_seqs = self._parse_fasta_all(efetch_plain.text)
            if plain_seqs:
                seq = self._sanitize_dna(plain_seqs[0])
                if seq:
                    return seq, f"web:ncbi:{nuccore_id}:fasta"

            raise ValueError(f"No usable sequence found on NCBI for '{gene_name}'.")

    def _fetch_by_ncbi_accession(self, accession: str) -> Tuple[str, str]:
        """
        Fetch a sequence directly from NCBI by accession number.
        Bypasses text search entirely — guarantees the correct record.
        """
        with httpx.Client(timeout=self._timeout_s, follow_redirects=True) as client:
            efetch = client.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params={
                    "db": "nuccore",
                    "id": accession,
                    "rettype": "fasta",
                    "retmode": "text",
                },
            )
            efetch.raise_for_status()
            seqs = self._parse_fasta_all(efetch.text)
            if seqs:
                seq = self._sanitize_dna(seqs[0])
                if seq:
                    return seq, f"web:ncbi:{accession}:direct"
        raise ValueError(f"No sequence returned for NCBI accession '{accession}'.")

    def _search_web_sequence(
        self, requested_name: str, element_type: str = "", is_backbone: bool = False
    ) -> Tuple[str, str]:
        # CDS / cargo: use dedicated clean-CDS fetcher.
        if element_type.lower() in (et.lower() for et in self._CDS_ELEMENT_TYPES):
            return self._search_web_cds(requested_name)

        # Build size-range filter to avoid returning whole plasmids for sub-elements,
        # or tiny fragments when a full backbone is expected.
        et_lower = element_type.lower()
        if is_backbone:
            size_filter = f"{self._BACKBONE_MIN_SIZE}:25000"
        else:
            size_filter = self._ELEMENT_SIZE_RANGES.get(et_lower, "")

        size_clause = f" AND {size_filter}[Sequence Length]" if size_filter else ""

        # Tailor keyword hints to element type so NCBI ranks the right records higher.
        if et_lower == "rep_origin":
            kw = "origin replication"
        elif et_lower in ("polya_signal", "terminator"):
            kw = "polyadenylation signal terminator"
        elif et_lower == "promoter":
            kw = "promoter"
        elif et_lower == "rbs":
            kw = "ribosome binding"
        else:
            kw = "promoter gene origin terminator marker"

        term = f"{requested_name}[Title]{size_clause} AND ({kw})"

        with httpx.Client(timeout=self._timeout_s, follow_redirects=True) as client:
            esearch = client.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params={"db": "nuccore", "retmode": "json", "retmax": 3, "term": term},
            )
            esearch.raise_for_status()
            ids = esearch.json().get("esearchresult", {}).get("idlist", [])

            # Retry without size filter if nothing found.
            if not ids and size_clause:
                term2 = f"{requested_name}[Title] AND ({kw})"
                esearch2 = client.get(
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                    params={"db": "nuccore", "retmode": "json", "retmax": 3, "term": term2},
                )
                esearch2.raise_for_status()
                ids = esearch2.json().get("esearchresult", {}).get("idlist", [])

            if not ids:
                raise ValueError(f"No NCBI hit for '{requested_name}'.")

            nuccore_id = ids[0]

            efetch = client.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params={"db": "nuccore", "id": nuccore_id, "rettype": "fasta", "retmode": "text"},
            )
            efetch.raise_for_status()
            seqs = self._parse_fasta_all(efetch.text)
            if not seqs:
                raise ValueError(f"No FASTA sequence for '{requested_name}'.")
            seq = self._sanitize_dna(seqs[0])
            if not seq:
                raise ValueError(f"Empty sequence for '{requested_name}'.")
            return seq, f"web:ncbi:{nuccore_id}"

    # ---------------------------------------------------------------------------
    # Mammalian backbone name canonicalisation
    # ---------------------------------------------------------------------------
    _MAMMALIAN_BACKBONE_ALIASES: Dict[str, str] = {
        "pcdna3.1": "pcDNA3.1",
        "pcdna3": "pcDNA3.1",
        "pcdna": "pcDNA3.1",
        "pcmv": "pCMV-Sport6",
        "pcmv-sport": "pCMV-Sport6",
        "pcag": "pCAGGS",
        "pcaggs": "pCAGGS",
        "paav": "pAAV-MCS",
        "paav-mcs": "pAAV-MCS",
        "plenti": "pLenti-CMV-GFP-Puro",
        "plncx": "pLNCX2",
        "pmscv": "pMSCV",
    }

    _MAMMALIAN_BACKBONE_HINTS = (
        "pcdna", "pcmv", "pcag", "pcagg", "paav", "plenti", "plncx", "pmscv",
    )

    @classmethod
    def _canonical_backbone_name(cls, name: str) -> str:
        """
        Map vague or partial mammalian backbone names to canonical identifiers
        known to Addgene/NCBI.  Returns the input unchanged when no match found.
        """
        if not name:
            return name
        n = name.strip().lower()
        # Direct alias lookup (exact or prefix match)
        for key, canonical in cls._MAMMALIAN_BACKBONE_ALIASES.items():
            if n == key or n.startswith(key):
                return canonical
        # Generic mammalian expression descriptor → default to pcDNA3.1
        mammalian_signals = ("mammalian", "mammal", "eukaryot", "cmv", "sv40")
        if sum(1 for kw in mammalian_signals if kw in n) >= 2:
            return "pcDNA3.1"
        return name

    # ---------------------------------------------------------------------------
    # Host detection
    # ---------------------------------------------------------------------------
    _MAMMALIAN_HOST_KEYWORDS = (
        "mammal", "mammalian", "hek", "hek293", "cho", "cos", "hela",
        "human", "mouse", "293", "vero", "jurkat", "eukaryot", "sf9", "insect",
    )

    @classmethod
    def _is_mammalian_host(cls, intent: IntentOutput) -> bool:
        """Return True when expression host is mammalian/eukaryotic."""
        for field in (intent.expression_host, intent.best_expression_host, intent.backbone):
            text = (field or "").lower()
            if any(kw in text for kw in cls._MAMMALIAN_HOST_KEYWORDS):
                return True
            if any(kw in text for kw in cls._MAMMALIAN_BACKBONE_HINTS):
                return True
        return False

    # ---------------------------------------------------------------------------
    # MCS selection
    # ---------------------------------------------------------------------------
    @classmethod
    def _default_mcs_for_backbone(cls, backbone_name: str) -> str:
        """
        Return the best-matching MCS name (from plasmid_MCS.csv) for a given
        backbone, so the correct polylinker is inserted when the backbone does
        not already contain one.
        """
        n = (backbone_name or "").lower()
        if any(k in n for k in ("pet", "prsf", "pacz", "prset", "pcold")):
            return "pET MCS"
        if any(k in n for k in ("pcdna", "pcmv", "pcag", "pef1", "plenti", "pmscv", "paav")):
            return "pcDNA3.1 MCS"
        if any(k in n for k in ("pbluescript", "pbs", "psk+", "pks+")):
            return "pBluescript MCS"
        # pUC19, pBR322, pGEM, pTZ, general bacterial → standard pUC19 polylinker
        return "pUC19 MCS"

    # ---------------------------------------------------------------------------
    # Feature detection
    # ---------------------------------------------------------------------------
    def _detect_features(self, sequence: str, backbone_name: str = "") -> Dict[str, bool]:
        """
        Scan backbone sequence *and* backbone name for common features.

        Returns:
            has_promoter            – any promoter signal detected
            has_mammalian_promoter  – CMV / SV40 / CAG / EF1 etc.
            has_polyA               – BGH / SV40 polyA signal detected
            has_resistance_marker   – antibiotic resistance gene detected
            has_mammalian_marker    – Neomycin/G418/Hygro/Puro detected
            has_cargo               – known fluorescent / reporter / target gene detected
            is_mammalian_backbone   – backbone name matches mammalian vector family
        """
        seq_lower = (sequence or "").lower()
        name_lower = (backbone_name or "").lower()
        combined = seq_lower + " " + name_lower

        mammalian_promoter_kw = ("cmv", "sv40", "cag", "ef1", "ef-1", "pgk", "ubc", "rous")
        bacterial_promoter_kw = ("plac", "ptac", "ptrc", "para", "ptet", "t7prom", "lacp", "arap")
        bacterial_seq_kw = ("taatacgactcactat",  # T7 core
                            "ttgacaattaatcatcg",  # tac
                            "ggaattgtgagcgctcac")  # lac
        marker_bacterial_kw = ("ampicillin", "bla", " kan ", "kanamycin",
                                "chloramphenicol", "spectinomycin", "tetracycline")
        marker_mammalian_kw = ("neomycin", "g418", "hygromycin", "puromycin",
                                "blasticidin", "zeocin")
        cargo_name_kw = ("gfp", "rfp", "mcherry", "egfp", "luciferase",
                          "tp53", "p53", "lacz", "beta-gal", "cargo", "insert")
        polya_kw = ("bgh", "sv40 polya", "sv40poly", "aataaa", "poly_a", "polya signal")

        has_mammalian_promoter = any(kw in combined for kw in mammalian_promoter_kw)
        has_bacterial_promoter = (
            any(kw in name_lower for kw in bacterial_promoter_kw)
            or any(kw in seq_lower for kw in bacterial_seq_kw)
        )
        has_polyA = any(kw in combined for kw in polya_kw)
        has_mammalian_marker = any(kw in combined for kw in marker_mammalian_kw)

        # Cargo: check name keywords + short seed match of gene name in sequence
        has_cargo_kw = any(kw in combined for kw in cargo_name_kw)

        # MCS detection: presence of ≥3 common restriction enzyme sites in the backbone
        # indicates a polylinker / multiple cloning site is already embedded.
        restriction_sites = (
            "gaattc",  # EcoRI
            "ggatcc",  # BamHI
            "aagctt",  # HindIII
            "ctgcag",  # PstI
            "gtcgac",  # SalI
            "gcatgc",  # SphI
            "cccggg",  # SmaI/XmaI
            "tctaga",  # XbaI
            "gagctc",  # SacI
            "ggtacc",  # KpnI
            "catatg",  # NdeI
            "ccatgg",  # NcoI
        )
        sites_found = sum(1 for site in restriction_sites if site in seq_lower)
        has_mcs = sites_found >= 3

        return {
            "has_promoter": has_mammalian_promoter or has_bacterial_promoter,
            "has_mammalian_promoter": has_mammalian_promoter,
            "has_polyA": has_polyA,
            "has_resistance_marker": (
                any(kw in combined for kw in marker_bacterial_kw)
                or has_mammalian_marker
            ),
            "has_mammalian_marker": has_mammalian_marker,
            "has_cargo": has_cargo_kw,
            "has_mcs": has_mcs,
            "is_mammalian_backbone": any(kw in name_lower for kw in self._MAMMALIAN_BACKBONE_HINTS),
            # Name lists (populated by _enrich_features_from_annotations when API data is available)
            "promoter_names": [],
            "terminator_names": [],
            "marker_names": [],
            "ori_names": [],
            "feature_source": "sequence_scan",
        }

    def _cargo_in_sequence(self, sequence: str, gene_symbol: str) -> bool:
        """Return True if gene_symbol (or its first 18 nt) appears in sequence."""
        if not sequence or not gene_symbol:
            return False
        seed = re.sub(r"[^ACGT]", "", gene_symbol.upper())[:18]
        if seed and seed in sequence.upper():
            return True
        return False

    # ---------------------------------------------------------------------------
    # Cargo removal
    # ---------------------------------------------------------------------------
    def _remove_cargo(
        self, sequence: str, gene_symbol: str = ""
    ) -> Tuple[str, List[str]]:
        """
        Remove cargo protein coding sequences from backbone.

        Strategy:
          1. If *gene_symbol* seed found in sequence, excise from ATG before it
             to nearest downstream stop codon.
          2. Fall back to the first generic ATG…stop ORF of 300-6000 bp.

        Returns:
            (cleaned_sequence, excision_log)   – log entries describe each
            excised region so callers can surface them as warnings.
        """
        if not sequence:
            return sequence, []

        excision_log: List[str] = []
        result = sequence

        # 1. Seed-based match for the requested gene (substring homology)
        if gene_symbol:
            seed = re.sub(r"[^ACGT]", "", gene_symbol.upper())[:18]
            if seed:
                idx = result.upper().find(seed)
                if idx >= 0:
                    # Walk back to the nearest ATG start codon (up to 300 bp upstream)
                    search_window = result[max(0, idx - 300): idx].upper()
                    atg_rel = search_window.rfind("ATG")
                    start = (idx - 300 + atg_rel) if atg_rel >= 0 else idx

                    # Find the nearest in-frame stop codon downstream
                    remainder = result[start:].upper()
                    stop_match = re.search(r"(?:TAA|TAG|TGA)", remainder[len(seed):])
                    end = start + len(seed) + (stop_match.end() if stop_match else 0)

                    excision_log.append(
                        f"gene_homology: excised [{start}:{end}] "
                        f"({end - start} bp) — seed matched '{gene_symbol}'"
                    )
                    result = result[:start] + result[end:]

        # 2. Generic CDS removal if still appears to have cargo
        cargo_name_kw = ("gfp", "rfp", "mcherry", "egfp", "luciferase", "tp53",
                          "p53", "lacz", "cargo")
        if any(kw in result.lower() for kw in cargo_name_kw):
            m = re.search(r"ATG[ACGT]{297,5997}(?:TAA|TAG|TGA)", result)
            if m:
                excision_log.append(
                    f"cds_pattern: excised [{m.start()}:{m.end()}] "
                    f"({m.end() - m.start()} bp) — generic ORF removal"
                )
                result = result[: m.start()] + result[m.end():]

        return result, excision_log

    def _resolve(
        self, slot: _Slot, is_backbone: bool = False
    ) -> Tuple[ResolvedElement, Optional[str]]:
        requested = slot.requested_name or f"default_{slot.name.lower().replace(' ', '_')}"

        # --- 1. Local dataset (returns sequence directly or an accession hint) ---
        from_dataset, acc_hint, acc_type = self._match_dataset(requested, slot.element_type)
        if from_dataset:
            return from_dataset.model_copy(update={"slot": slot.name}), None

        seq: str = ""
        source: str = ""

        # --- 2. Direct NCBI accession fetch (from CSV metadata row) ---
        if acc_hint and acc_type == "ncbi":
            try:
                seq, source = self._fetch_by_ncbi_accession(acc_hint)
                self._logger.info(
                    "[CONSTRUCTION] fetched '%s' via accession %s (%d bp)",
                    requested, acc_hint, len(seq),
                )
            except Exception as acc_exc:
                self._logger.info(
                    "[CONSTRUCTION] accession fetch failed for %s (%s): %s",
                    requested, acc_hint, acc_exc,
                )

        # Direct Addgene fetch when accession points there and NCBI fetch missed.
        if not seq and acc_hint and acc_type == "addgene":
            try:
                seq, source = self._search_addgene_sequence(acc_hint)
            except Exception as ag_exc:
                self._logger.info(
                    "[CONSTRUCTION] addgene direct fetch failed for %s (%s): %s",
                    requested, acc_hint, ag_exc,
                )

        # --- 3. Addgene text-search fallback (by name) ---
        if not seq:
            try:
                seq, source = self._search_addgene_sequence(requested)
            except Exception as addgene_exc:
                self._logger.info(
                    "[CONSTRUCTION] addgene fallback failed for %s: %s", requested, addgene_exc
                )

        # --- 4. NCBI text-search fallback (element-type-aware, size-filtered) ---
        if not seq:
            try:
                seq, source = self._search_web_sequence(
                    requested, slot.element_type, is_backbone=is_backbone
                )
            except Exception as web_exc:
                self._logger.info(
                    "[CONSTRUCTION] web fallback failed for %s: %s", requested, web_exc
                )

        # --- 5. Last resort placeholder (disabled in strict mode) ---
        if not seq:
            if self._strict:
                raise ValueError(
                    f"[strict mode] No sequence found for '{requested}' "
                    f"(type={slot.element_type}). "
                    "Check dataset CSV or disable strict mode to allow placeholders."
                )
            seq = self._placeholder_seq(requested, slot.element_type)
            source = "web:fallback-placeholder"

        web_hit = ResolvedElement(
            slot=slot.name,
            requested_name=requested,
            element_type=slot.element_type,
            sequence=seq,
            source=source,
            from_web=source != "web:fallback-placeholder",
        )
        return web_hit, WEB_FALLBACK_WARNING

    def _template_slots(
        self,
        intent: IntentOutput,
        features: Optional[Dict] = None,
        is_mammalian: bool = False,
        backbone_name: str = "",
        skip_reasons: Optional[List[str]] = None,
    ) -> List[_Slot]:
        """
        Build ordered slot list from construct type, backbone features, and host.

        Deduplication rules
        -------------------
        * Backbone already has a promoter   → skip Promoter slot(s).
        * Backbone already has a polyA/terminator → skip Terminator/polyA slot.
        * Backbone already has a resistance marker → skip Selection Marker slot.
        * Mammalian host                    → omit all RBS slots; default to
          Neomycin/G418 for the selection marker.

        skip_reasons (optional mutable list): populated with human-readable
        descriptions of every slot omission so callers can surface them as
        frontend warnings / SSE events.
        """
        f = features or {}
        skip_promoter = f.get("has_promoter", False)
        skip_polya = f.get("has_polyA", False)
        skip_marker = f.get("has_resistance_marker", False)
        skip_mcs = f.get("has_mcs", False)   # backbone already has a polylinker
        src = f.get("feature_source", "sequence_scan")

        # Human-readable name lists (populated by Addgene/GenBank annotation)
        promoter_names_str = ", ".join(f.get("promoter_names") or [])
        terminator_names_str = ", ".join(f.get("terminator_names") or [])
        marker_names_str = ", ".join(f.get("marker_names") or [])

        def _log_skip(slot_label: str, reason: str) -> None:
            msg = f"Skipping {slot_label}: {reason} (detected via {src})"
            self._logger.info("[CONSTRUCTION] %s", msg)
            if skip_reasons is not None:
                skip_reasons.append(msg)

        # Mammalian default marker
        def _marker_name() -> str:
            if intent.selection_marker:
                return intent.selection_marker
            return "Neomycin (G418)" if is_mammalian else None

        # Terminator / polyA slot label
        def _term_slot(label: str) -> _Slot:
            name = intent.polyA or intent.terminator
            slot_type = "polyA_signal" if is_mammalian else "terminator"
            return _Slot(label, slot_type, name or ("BGH polyA" if is_mammalian else None))

        # -----------------------------------------------------------------------
        # Fusion construct
        # -----------------------------------------------------------------------
        if intent.construct_type == "fusion":
            slots: List[_Slot] = [
                _Slot("Origin", "rep_origin", intent.origin_of_replication),
            ]
            if not skip_promoter:
                slots.append(_Slot("Promoter1", "promoter", intent.promoter))
            if not is_mammalian:
                slots.append(_Slot("RBS", "RBS", intent.rbs or "RBS_B0034"))
            slots.extend([
                _Slot("Gene A", "CDS", intent.gene_a or intent.gene_symbol),
                _Slot("Linker Sequence", "misc_feature", intent.linker_sequence or "GGGGSx3"),
                _Slot("Gene B", "CDS", intent.gene_b or "fusion_partner"),
            ])
            if not skip_polya:
                slots.append(_term_slot("Terminator / polyA"))
            if not skip_marker:
                slots.append(_Slot("Selection Marker", "CDS", _marker_name()))
            if intent.regulatory_element:
                slots.append(_Slot("Regulatory Element", "misc_feature", intent.regulatory_element))
            return slots

        # -----------------------------------------------------------------------
        # Multi-cassette fusion
        # -----------------------------------------------------------------------
        if intent.construct_type == "multi_cassette_fusion":
            slots = [
                _Slot("Origin", "rep_origin", intent.origin_of_replication),
            ]
            if not skip_promoter:
                slots.append(_Slot("Promoter1", "promoter", intent.promoter))
            if not is_mammalian:
                slots.append(_Slot("RBS1", "RBS", intent.rbs or "RBS_B0034"))
            slots.extend([
                _Slot("Gene A", "CDS", intent.gene_a or intent.gene_symbol),
                _Slot("Stop", "misc_feature", "stop_codon"),
                _Slot("Promoter2", "promoter", intent.promoter2 or intent.promoter),
            ])
            if not is_mammalian:
                slots.append(_Slot("RBS2", "RBS", intent.rbs2 or intent.rbs or "RBS_B0032"))
            slots.append(_Slot("Gene B", "CDS", intent.gene_b or "second_gene"))
            if not skip_polya:
                slots.append(_term_slot("Terminator / polyA"))
            if not skip_marker:
                slots.append(_Slot("Selection Marker", "CDS", _marker_name()))
            return slots

        # -----------------------------------------------------------------------
        # Synthetic circuit
        # -----------------------------------------------------------------------
        if intent.construct_type == "synthetic_circuit":
            slots = [
                _Slot("Origin of Replication", "rep_origin", intent.origin_of_replication),
            ]
            # If explicit components list provided by intent agent, use them directly.
            if intent.circuit_components:
                for comp in intent.circuit_components:
                    slots.append(
                        _Slot(comp.name, comp.element_type, comp.sequence_name or comp.name)
                    )
            else:
                # Fallback: build a minimal circuit from the high-level fields.
                # Pattern: Promoter → RBS → Regulator(s) → Terminator
                #          → [second cassette: sensor promoter → RBS → reporter → Terminator]
                if not skip_promoter:
                    slots.append(_Slot("Promoter", "promoter", intent.promoter))
                if not is_mammalian:
                    slots.append(_Slot("RBS", "RBS", intent.rbs or "RBS_B0034"))
                # Regulatory genes (e.g. TetR, LacI)
                for reg in intent.circuit_regulators:
                    slots.append(_Slot(f"Regulator ({reg})", "CDS", reg))
                if not skip_polya:
                    slots.append(_term_slot("Terminator 1"))
                # Reporter cassette
                if intent.output_reporters:
                    for i, reporter in enumerate(intent.output_reporters):
                        tag = f" {i + 1}" if len(intent.output_reporters) > 1 else ""
                        if not skip_promoter:
                            slots.append(
                                _Slot(f"Sensor Promoter{tag}", "promoter",
                                      intent.promoter2 or intent.promoter)
                            )
                        if not is_mammalian:
                            slots.append(_Slot(f"RBS{tag}", "RBS", intent.rbs2 or intent.rbs or "RBS_B0034"))
                        slots.append(_Slot(f"Reporter{tag} ({reporter})", "CDS", reporter))
                        if not skip_polya:
                            slots.append(_term_slot(f"Terminator{tag}"))
                elif intent.gene_symbol:
                    # No explicit reporters — treat gene_symbol as the output
                    slots.append(_Slot("Gene of Interest", "CDS", intent.gene_symbol))
                    if not skip_polya:
                        slots.append(_term_slot("Terminator"))
            if not skip_marker:
                slots.append(_Slot("Selection Marker", "CDS", _marker_name()))
            if not skip_mcs:
                mcs_name = intent.cloning_site or self._default_mcs_for_backbone(backbone_name)
                slots.append(_Slot("MCS/Cloning Site", "misc_feature", mcs_name))
            return slots

        # -----------------------------------------------------------------------
        # Default single-gene cassette
        # -----------------------------------------------------------------------
        slots = [
            _Slot("Origin of Replication", "rep_origin", intent.origin_of_replication),
        ]
        if not skip_promoter:
            slots.append(_Slot("Promoter", "promoter", intent.promoter))
        if not is_mammalian:
            slots.append(_Slot("RBS", "RBS", intent.rbs or "RBS_B0034"))
        slots.append(_Slot("Gene of Interest", "CDS", intent.gene_symbol))
        if not skip_polya:
            slots.append(_term_slot("Terminator / polyA"))
        if not skip_marker:
            slots.append(_Slot("Selection Marker Cassette", "CDS", _marker_name()))
        if intent.regulatory_element:
            slots.append(_Slot("Regulatory Element", "misc_feature", intent.regulatory_element))
        # MCS: skip if backbone already has a polylinker (≥3 restriction sites detected).
        # Otherwise always add one — use the user-specified site or derive from backbone type.
        if not skip_mcs:
            mcs_name = intent.cloning_site or self._default_mcs_for_backbone(backbone_name)
            slots.append(_Slot("MCS/Cloning Site", "misc_feature", mcs_name))
        return slots

    def run(
        self,
        inp: PlasmidConstructionInput,
        progress_cb: Optional[Callable[[Dict], None]] = None,
    ) -> PlasmidConstructionOutput:
        """
        progress_cb, if provided, is called with a dict after each significant step:
            {"step": "backbone"|"element"|"done", "label": str, "bp": int, "source": str}
        """
        def _emit(event: Dict) -> None:
            if progress_cb:
                try:
                    progress_cb(event)
                except Exception:
                    pass

        self._logger.info("[CONSTRUCTION] INPUT %s", inp.model_dump())
        warnings: List[str] = []
        elements: List[ResolvedElement] = []

        # Step 1: Canonicalise backbone name, detect host type
        is_mammalian = self._is_mammalian_host(inp.intent)
        raw_backbone_name = inp.intent.backbone or ("pcDNA3.1" if is_mammalian else "pUC19")
        backbone_name = self._canonical_backbone_name(raw_backbone_name)
        if backbone_name != raw_backbone_name:
            self._logger.info(
                "[CONSTRUCTION] Backbone name canonicalised: '%s' → '%s'",
                raw_backbone_name, backbone_name,
            )
            warnings.append(
                f"Backbone name normalised from '{raw_backbone_name}' to '{backbone_name}'."
            )

        self._logger.info(
            "[CONSTRUCTION] Resolving backbone '%s' (mammalian=%s)", backbone_name, is_mammalian
        )
        _emit({"step": "resolving", "label": "Backbone", "name": backbone_name})
        backbone_slot = _Slot("Backbone", "rep_origin", backbone_name)
        backbone, backbone_warn = self._resolve(backbone_slot, is_backbone=True)
        if backbone_warn:
            warnings.append(backbone_warn)

        self._logger.info("[CONSTRUCTION] Backbone sequence length: %d bp", len(backbone.sequence))
        _emit({
            "step": "backbone",
            "label": "Backbone",
            "name": backbone_name,
            "bp": len(backbone.sequence),
            "source": backbone.source,
        })

        # Step 2: Detect features in backbone (sequence + name)
        features = self._detect_features(backbone.sequence, backbone_name)
        self._logger.info(
            "[CONSTRUCTION] Backbone features: mammalian_promoter=%s, polyA=%s, "
            "mammalian_marker=%s, cargo=%s",
            features.get("has_mammalian_promoter"),
            features.get("has_polyA"),
            features.get("has_mammalian_marker"),
            features.get("has_cargo"),
        )

        # Step 3: If cargo detected, excise it before inserting new gene
        cleaned_backbone_seq = backbone.sequence
        cargo_detected = features.get("has_cargo") or self._cargo_in_sequence(
            backbone.sequence, inp.intent.gene_symbol
        )
        if cargo_detected:
            self._logger.info("[CONSTRUCTION] Cargo detected in backbone; removing...")
            cleaned_backbone_seq, excision_log = self._remove_cargo(
                backbone.sequence, inp.intent.gene_symbol
            )
            self._logger.info(
                "[CONSTRUCTION] Cleaned backbone: %d bp → %d bp",
                len(backbone.sequence), len(cleaned_backbone_seq),
            )
            for entry in excision_log:
                msg = f"Cargo excised from backbone — {entry}"
                self._logger.info("[CONSTRUCTION] %s", msg)
                warnings.append(msg)

        # Create cleaned backbone element
        backbone_cleaned = backbone.model_copy(update={"sequence": cleaned_backbone_seq})
        elements.append(backbone_cleaned)

        # Step 4: Build template slots (deduplication + organism-aware)
        slots = self._template_slots(
            inp.intent, features=features, is_mammalian=is_mammalian, backbone_name=backbone_name
        )
        self._logger.info("[CONSTRUCTION] Template slots: %d", len(slots))
        
        # Step 5: Resolve remaining elements
        for slot in slots:
            _emit({"step": "resolving", "label": slot.name, "name": slot.requested_name or slot.name})
            resolved, warn = self._resolve(slot)
            elements.append(resolved)
            if warn and warn not in warnings:
                warnings.append(warn)
            _emit({
                "step": "element",
                "label": slot.name,
                "name": resolved.requested_name,
                "bp": len(resolved.sequence),
                "source": resolved.source,
                "from_web": resolved.from_web,
            })

        # Step 6: Assemble final construct (backbone + other elements)
        seq = cleaned_backbone_seq + "".join(e.sequence for e in elements[1:])
        
        out = PlasmidConstructionOutput(
            template=inp.intent.construct_type,
            elements=elements,
            construct_sequence=seq,
            frontend_warnings=warnings,
        )
        self._logger.info("[CONSTRUCTION] OUTPUT %s", out.model_dump())
        return out