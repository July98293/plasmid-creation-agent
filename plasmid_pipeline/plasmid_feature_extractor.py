from __future__ import annotations

import io
import re
from dataclasses import dataclass
from typing import Iterable, Optional

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


@dataclass(frozen=True)
class ExtractedFeature:
    name: str
    feature_type: str
    sequence: str
    start_1: int
    end_1: int
    strand: int  # 1 or -1 (or 0 if unknown)


class PlasmidFeatureExtractor:
    """
    Extracts a single annotated feature from a GenBank record.

    This is intentionally strict about extracting ONLY the requested region,
    never the whole plasmid sequence.
    """

    def _parse_first_record(self, genbank_text: str):
        handle = io.StringIO(genbank_text)
        return next(SeqIO.parse(handle, "genbank"))

    def _qualifier_text(self, feat: SeqFeature) -> str:
        parts: list[str] = [feat.type]
        q = feat.qualifiers or {}
        for k, v in q.items():
            if isinstance(v, list):
                parts.extend(str(x) for x in v)
            else:
                parts.append(str(v))
        return " ".join(parts).lower()

    def _matches_any(self, haystack: str, needles: Iterable[str]) -> bool:
        h = haystack.lower()
        for n in needles:
            n = (n or "").strip().lower()
            if not n:
                continue
            if n in h:
                return True
        return False

    def extract_by_keywords(
        self,
        *,
        genbank_text: str,
        desired_feature_types: set[str],
        keywords: list[str],
        prefer_exact_word: bool = False,
        name_for_output: Optional[str] = None,
    ) -> Optional[ExtractedFeature]:
        """
        Find the best matching feature of one of `desired_feature_types` whose
        qualifiers mention any of `keywords`, then return the extracted sequence.
        """
        record = self._parse_first_record(genbank_text)

        candidates: list[SeqFeature] = [
            f for f in getattr(record, "features", []) if getattr(f, "type", None) in desired_feature_types
        ]
        if not candidates:
            return None

        # First pass: match by keyword substring.
        matched: list[SeqFeature] = []
        for feat in candidates:
            text = self._qualifier_text(feat)
            if prefer_exact_word:
                for kw in keywords:
                    kw2 = (kw or "").strip().lower()
                    if not kw2:
                        continue
                    if re.search(rf"\b{re.escape(kw2)}\b", text):
                        matched.append(feat)
                        break
            else:
                if self._matches_any(text, keywords):
                    matched.append(feat)

        if not matched:
            return None

        # Prefer the shortest match (most likely the part, not a whole cassette).
        def feat_len(f: SeqFeature) -> int:
            try:
                return int(len(f.location))
            except Exception:
                return 10**9

        best = sorted(matched, key=feat_len)[0]

        try:
            seq = str(best.extract(record.seq)).upper()
        except Exception:
            return None

        # Coordinates: BioPython uses 0-based, end-exclusive.
        parts = list(best.location.parts) if hasattr(best.location, "parts") else [best.location]
        starts = [int(p.start) for p in parts]
        ends = [int(p.end) for p in parts]
        start_1 = min(starts) + 1
        end_1 = max(ends)

        strand = int(getattr(best.location, "strand", 0) or 0)
        if strand not in (-1, 0, 1):
            strand = 0

        return ExtractedFeature(
            name=name_for_output or (keywords[0] if keywords else "extracted_feature"),
            feature_type=best.type,
            sequence=seq,
            start_1=start_1,
            end_1=end_1,
            strand=strand,
        )

