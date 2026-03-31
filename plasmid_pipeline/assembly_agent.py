from __future__ import annotations

import io
import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from logging_utils import get_conversation_logger
from .models import (
    Annotation,
    AssemblyFragment,
    AssemblyInput,
    AssemblyOutput,
    WarningMessage,
)

_RC_TABLE = str.maketrans("ATGCN", "TACGN")
logger = logging.getLogger(__name__)


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


# Anything in these classes is considered part of the expression cassette logic.
_CASSETTE_FEATURE_TYPES = {
    "promoter",
    "terminator",
    "polyA_signal",
    "polyA signal",
    "sig_peptide",
    "5'UTR",
    "3'UTR",
    "misc_feature",
    "CDS",
    "kozak",
    "tag",
}


@dataclass
class BackboneFeatureHit:
    name: str
    type: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    sequence: str
    source: str  # "genbank" or "sequence_scan"


class AssemblyAgent:
    """
    Improved biologically-aware assembly planner with deduplication logic.

    Main changes vs previous version:
    - Reads backbone GenBank annotations if present.
    - Also scans the raw backbone sequence directly for cassette features.
    - If a promoter/CDS/terminator/etc. already exists in the backbone and is
      also present in the insert, removes the old one and inserts the new one.
    - Uses sequence-level conflict detection instead of relying only on GenBank.
    - Filters cassette features to avoid re-inserting features already in backbone.
    - Rebuilds cassette with only non-duplicate features.
    - Helps downstream "feature not found" situations by recognizing features
      already present in the backbone even when they were not explicitly resolved
      upstream.
    """

    def __init__(self) -> None:
        self._logger = get_conversation_logger()

    # ------------------------------------------------------------------
    # Basic helpers
    # ------------------------------------------------------------------

    def _clean_seq(self, seq: str) -> str:
        """Remove non-standard nucleotides from sequence."""
        return "".join(ch for ch in (seq or "").upper() if ch in {"A", "T", "G", "C", "N"})

    def _classify_host(self, expression_host: str) -> str:
        """Classify expression host as mammalian, bacterial, or unknown."""
        host = (expression_host or "").lower()
        for kw in ("hek", "cho", "cos", "vero", "hela", "293", "mammalian", "human", "mouse", "rat"):
            if kw in host:
                return "mammalian"
        for kw in ("e.coli", "ecoli", "bl21", "dh5", "bacteria", "bacterial", "prokaryot"):
            if kw in host:
                return "bacterial"
        return "unknown"

    def _pick_method(self, pref: str) -> str:
        """Normalize assembly method preference to canonical name."""
        p = (pref or "").lower()
        if "gibson" in p:
            return "Gibson"
        if "golden" in p:
            return "Golden Gate"
        if "restrict" in p or "digest" in p:
            return "Restriction"
        if "synth" in p:
            return "Synthesis"
        return "Gibson"

    def _normalize_feature_type(self, feat_type: str) -> str:
        """Normalize feature type names for consistency."""
        t = (feat_type or "").strip()
        if t.lower() == "kozak":
            return "kozak"
        if t.lower() in {"polya_signal", "polya signal"}:
            return "polyA_signal"
        return t

    # ------------------------------------------------------------------
    # GenBank parsing
    # ------------------------------------------------------------------

    def _parse_genbank_features(self, genbank_text: str) -> List[BackboneFeatureHit]:
        """
        Parse raw GenBank text and return feature hits from existing annotations.
        """
        hits: List[BackboneFeatureHit] = []
        if not genbank_text:
            return hits

        try:
            record: SeqRecord = SeqIO.read(io.StringIO(genbank_text), "genbank")
        except Exception as exc:
            self._logger.warning("[ASSEMBLY] GenBank parse failed: %s", exc)
            return hits

        full_seq = str(record.seq).upper()

        for feat in record.features:
            if feat.type == "source":
                continue

            name = (
                feat.qualifiers.get("label", [""])[0]
                or feat.qualifiers.get("gene", [""])[0]
                or feat.qualifiers.get("product", [""])[0]
                or feat.type
            )

            try:
                start = int(feat.location.start)
                end = int(feat.location.end)
                seq = full_seq[start:end]
            except Exception:
                continue

            hits.append(
                BackboneFeatureHit(
                    name=name,
                    type=self._normalize_feature_type(feat.type),
                    start=start,
                    end=end,
                    sequence=seq,
                    source="genbank",
                )
            )

        self._logger.info(
            "[ASSEMBLY] Parsed GenBank backbone features: %s",
            [f"{h.name}:{h.type}[{h.start}:{h.end}]" for h in hits],
        )
        return hits

    # ------------------------------------------------------------------
    # Cassette building
    # ------------------------------------------------------------------

    def _build_cassette(
        self,
        features,
        cds_sequence: str,
        gene_symbol: str,
        warnings: List[WarningMessage],
    ) -> Tuple[str, str, str, List[str], List[Annotation], List[Dict]]:
        """
        Build the expression cassette:
            promoter + kozak + [n_tags] + cds + [unspecified_tags] + [c_tags] + terminator

        Returns:
            cassette,
            clean_cds,
            promoter_seq,
            feature_order,
            cassette_annotations,
            cassette_feature_dicts
        """
        promoters = [f for f in features if f.type == "promoter"]
        kozak_feats = [f for f in features if f.type == "kozak"]
        tag_feats = [f for f in features if f.type == "tag"]
        terminators = [f for f in features if f.type == "terminator"]

        n_tags = [t for t in tag_feats if (t.position_hint or "").upper() == "N"]
        c_tags = [t for t in tag_feats if (t.position_hint or "").upper() == "C"]
        unspecified_tags = [t for t in tag_feats if (t.position_hint or "").upper() not in ("N", "C")]

        if unspecified_tags:
            warnings.append(
                WarningMessage(
                    code="TAG_POSITION_UNSPECIFIED",
                    message=(
                        "One or more tag features have no position_hint. "
                        "Unspecified tags will be placed after the CDS by default."
                    ),
                )
            )

        parts: List[str] = []
        feature_order: List[str] = []
        annotations: List[Annotation] = []
        cassette_feature_dicts: List[Dict] = []
        cursor = 1

        def _append(name: str, seq: str, feat_type: str, label: str) -> None:
            nonlocal cursor
            clean = self._clean_seq(seq)
            if not clean:
                return

            parts.append(clean)
            annotations.append(
                Annotation(name=name, start=cursor, end=cursor + len(clean) - 1, type=feat_type)
            )
            cassette_feature_dicts.append(
                {
                    "name": name,
                    "type": self._normalize_feature_type(feat_type),
                    "start": cursor - 1,                # cassette-local 0-based
                    "end": cursor - 1 + len(clean),     # cassette-local 0-based exclusive
                    "sequence": clean,
                }
            )
            feature_order.append(label)
            cursor += len(clean)

        # Promoter
        if promoters:
            p = promoters[0]
            _append(p.name, p.sequence, "promoter", f"{p.name}_promoter")
            self._logger.info("[ASSEMBLY] promoter (%s): %d bp", p.name, len(self._clean_seq(p.sequence)))
        else:
            warnings.append(
                WarningMessage(
                    code="NO_PROMOTER",
                    message="No promoter feature resolved; construct may not be expressed.",
                )
            )
            self._logger.info("[ASSEMBLY] promoter: 0 bp (not resolved)")

        # Kozak
        if kozak_feats:
            k = kozak_feats[0]
            _append("Kozak", k.sequence, "kozak", "Kozak")
            self._logger.info("[ASSEMBLY] kozak: %d bp", len(self._clean_seq(k.sequence)))
        else:
            self._logger.info("[ASSEMBLY] kozak: 0 bp")

        # N-terminal tags
        for tag in n_tags:
            _append(f"{tag.name}_tag", tag.sequence, "tag", f"{tag.name}_N_tag")

        # CDS
        clean_cds = self._clean_seq(cds_sequence)
        if not clean_cds:
            warnings.append(
                WarningMessage(code="NO_CDS", message="No CDS sequence available; construct will be incomplete.")
            )
        else:
            if len(clean_cds) % 3 != 0:
                warnings.append(
                    WarningMessage(code="CDS_NOT_MULTIPLE_OF_3", message="CDS length is not divisible by 3.")
                )
            if not clean_cds.startswith("ATG"):
                warnings.append(
                    WarningMessage(code="CDS_NO_START_CODON", message="CDS does not start with ATG.")
                )
            if clean_cds[-3:] not in {"TAA", "TAG", "TGA"}:
                warnings.append(
                    WarningMessage(code="CDS_NO_STOP_CODON", message="CDS does not end with a standard stop codon.")
                )
            _append(f"{gene_symbol}_CDS", clean_cds, "CDS", f"{gene_symbol}_CDS")
            self._logger.info("[ASSEMBLY] cds: %d bp", len(clean_cds))

        # Unspecified tags
        for tag in unspecified_tags:
            _append(f"{tag.name}_tag", tag.sequence, "tag", f"{tag.name}_tag")

        # C-terminal tags
        for tag in c_tags:
            _append(f"{tag.name}_tag", tag.sequence, "tag", f"{tag.name}_C_tag")

        # Terminator
        if terminators:
            t = terminators[0]
            _append(t.name, t.sequence, "terminator", t.name)
            self._logger.info("[ASSEMBLY] terminator (%s): %d bp", t.name, len(self._clean_seq(t.sequence)))
        else:
            self._logger.info("[ASSEMBLY] terminator: 0 bp (not resolved)")

        cassette = "".join(parts)
        self._logger.info("[ASSEMBLY] cassette total: %d bp", len(cassette))

        promoter_seq = self._clean_seq(promoters[0].sequence) if promoters else ""

        # Structural validation
        if clean_cds and clean_cds not in cassette:
            raise ValueError("Cassette missing CDS")
        if promoter_seq and promoter_seq not in cassette:
            raise ValueError("Cassette missing promoter")
        if promoter_seq and len(cassette) < len(clean_cds) + 200:
            warnings.append(
                WarningMessage(
                    code="CASSETTE_TOO_SHORT",
                    message=(
                        f"Cassette is {len(cassette)} bp but expected at least "
                        f"{len(clean_cds) + 200} bp with flanking elements. "
                        "Promoter or terminator sequences may be incomplete."
                    ),
                )
            )

        return cassette, clean_cds, promoter_seq, feature_order, annotations, cassette_feature_dicts

    # ------------------------------------------------------------------
    # Feature filtering (deduplication)
    # ------------------------------------------------------------------

    def _filter_cassette_features(
        self,
        cassette_feature_dicts: List[Dict],
        backbone: str,
        backbone_hits: List[BackboneFeatureHit],
    ) -> Tuple[List[Dict], List[WarningMessage]]:
        """
        Remove features from cassette that are already present in backbone.

        Compare each cassette feature (type + sequence) against backbone features.
        If a match is found, the feature is skipped and a warning is added.

        Returns:
            filtered_features, skip_warnings
        """
        warnings: List[WarningMessage] = []
        filtered: List[Dict] = []

        # Create a set of (type, sequence) tuples already in backbone
        backbone_feature_seqs: Dict[Tuple[str, str], BackboneFeatureHit] = {}
        for h in backbone_hits:
            key = (self._normalize_feature_type(h.type), self._clean_seq(h.sequence))
            # Keep the sequence_scan hit if available (more reliable than genbank)
            if key not in backbone_feature_seqs or h.source == "sequence_scan":
                backbone_feature_seqs[key] = h

        self._logger.info(
            "[ASSEMBLY] Backbone feature fingerprints: %d unique features",
            len(backbone_feature_seqs),
        )

        for feat in cassette_feature_dicts:
            feat_type = self._normalize_feature_type(feat["type"])
            feat_seq = self._clean_seq(feat["sequence"])

            # Check if this feature already exists in backbone
            key = (feat_type, feat_seq)
            if key in backbone_feature_seqs:
                backbone_hit = backbone_feature_seqs[key]
                warnings.append(
                    WarningMessage(
                        code="FEATURE_ALREADY_IN_BACKBONE",
                        message=(
                            f"Feature '{feat['name']}' ({feat_type}, {len(feat_seq)} bp) is already present "
                            f"in the backbone at position [{backbone_hit.start}:{backbone_hit.end}] "
                            f"(source: {backbone_hit.source}) and will not be re-inserted."
                        ),
                    )
                )
                self._logger.info(
                    "[ASSEMBLY] Skipping cassette feature '%s' (%s, %d bp) — already in backbone "
                    "at [%d:%d] from %s",
                    feat["name"],
                    feat_type,
                    len(feat_seq),
                    backbone_hit.start,
                    backbone_hit.end,
                    backbone_hit.source,
                )
                continue

            filtered.append(feat)

        self._logger.info(
            "[ASSEMBLY] Filtered cassette features: %d → %d (removed %d duplicates)",
            len(cassette_feature_dicts),
            len(filtered),
            len(cassette_feature_dicts) - len(filtered),
        )

        return filtered, warnings

    def _rebuild_cassette_from_filtered(
        self,
        features,
        filtered_cassette_feature_dicts: List[Dict],
        cds_sequence: str,
        gene_symbol: str,
        warnings: List[WarningMessage],
    ) -> Tuple[str, str, str, List[str], List[Annotation], List[Dict]]:
        """
        Rebuild the expression cassette using only the filtered features.

        Maintains the same structural order:
            promoter + kozak + [n_tags] + cds + [unspecified_tags] + [c_tags] + terminator

        Any feature that was filtered out (already in backbone) is not included.

        Returns:
            cassette, clean_cds, promoter_seq, feature_order, annotations, cassette_feature_dicts
        """

        # Create a lookup set for filtered feature names
        filtered_names: set[str] = {f["name"] for f in filtered_cassette_feature_dicts}

        # Reconstruct feature lists from filtered data only
        promoters = [f for f in features if f.type == "promoter" and f.name in filtered_names]
        kozak_feats = [f for f in features if f.type == "kozak" and f.name in filtered_names]
        tag_feats = [f for f in features if f.type == "tag" and f.name in filtered_names]
        terminators = [f for f in features if f.type == "terminator" and f.name in filtered_names]

        n_tags = [t for t in tag_feats if (t.position_hint or "").upper() == "N"]
        c_tags = [t for t in tag_feats if (t.position_hint or "").upper() == "C"]
        unspecified_tags = [t for t in tag_feats if (t.position_hint or "").upper() not in ("N", "C")]

        if unspecified_tags:
            warnings.append(
                WarningMessage(
                    code="TAG_POSITION_UNSPECIFIED",
                    message=(
                        "One or more tag features have no position_hint. "
                        "Unspecified tags will be placed after the CDS by default."
                    ),
                )
            )

        parts: List[str] = []
        feature_order: List[str] = []
        annotations: List[Annotation] = []
        cassette_feature_dicts_rebuilt: List[Dict] = []
        cursor = 1

        def _append(name: str, seq: str, feat_type: str, label: str) -> None:
            nonlocal cursor
            clean = self._clean_seq(seq)
            if not clean:
                return

            parts.append(clean)
            annotations.append(
                Annotation(name=name, start=cursor, end=cursor + len(clean) - 1, type=feat_type)
            )
            cassette_feature_dicts_rebuilt.append(
                {
                    "name": name,
                    "type": self._normalize_feature_type(feat_type),
                    "start": cursor - 1,
                    "end": cursor - 1 + len(clean),
                    "sequence": clean,
                }
            )
            feature_order.append(label)
            cursor += len(clean)

        # Promoter
        if promoters:
            p = promoters[0]
            _append(p.name, p.sequence, "promoter", f"{p.name}_promoter")
            self._logger.info("[ASSEMBLY] promoter (%s): %d bp", p.name, len(self._clean_seq(p.sequence)))
        else:
            self._logger.info("[ASSEMBLY] promoter: 0 bp (not in cassette after filtering)")

        # Kozak
        if kozak_feats:
            k = kozak_feats[0]
            _append("Kozak", k.sequence, "kozak", "Kozak")
            self._logger.info("[ASSEMBLY] kozak: %d bp", len(self._clean_seq(k.sequence)))

        # N-terminal tags
        for tag in n_tags:
            _append(f"{tag.name}_tag", tag.sequence, "tag", f"{tag.name}_N_tag")

        # CDS
        clean_cds = self._clean_seq(cds_sequence)
        if not clean_cds:
            warnings.append(
                WarningMessage(code="NO_CDS", message="No CDS sequence available; construct will be incomplete.")
            )
        else:
            if len(clean_cds) % 3 != 0:
                warnings.append(
                    WarningMessage(code="CDS_NOT_MULTIPLE_OF_3", message="CDS length is not divisible by 3.")
                )
            if not clean_cds.startswith("ATG"):
                warnings.append(
                    WarningMessage(code="CDS_NO_START_CODON", message="CDS does not start with ATG.")
                )
            if clean_cds[-3:] not in {"TAA", "TAG", "TGA"}:
                warnings.append(
                    WarningMessage(code="CDS_NO_STOP_CODON", message="CDS does not end with a standard stop codon.")
                )
            _append(f"{gene_symbol}_CDS", clean_cds, "CDS", f"{gene_symbol}_CDS")
            self._logger.info("[ASSEMBLY] cds: %d bp", len(clean_cds))

        # Unspecified tags
        for tag in unspecified_tags:
            _append(f"{tag.name}_tag", tag.sequence, "tag", f"{tag.name}_tag")

        # C-terminal tags
        for tag in c_tags:
            _append(f"{tag.name}_tag", tag.sequence, "tag", f"{tag.name}_C_tag")

        # Terminator
        if terminators:
            t = terminators[0]
            _append(t.name, t.sequence, "terminator", t.name)
            self._logger.info("[ASSEMBLY] terminator (%s): %d bp", t.name, len(self._clean_seq(t.sequence)))

        cassette = "".join(parts)
        self._logger.info("[ASSEMBLY] cassette total (after filtering): %d bp", len(cassette))

        promoter_seq = self._clean_seq(promoters[0].sequence) if promoters else ""

        # Structural validation
        if clean_cds and clean_cds not in cassette:
            raise ValueError("Cassette missing CDS after filtering")
        if promoter_seq and promoter_seq not in cassette:
            raise ValueError("Cassette missing promoter after filtering")
        if promoter_seq and len(cassette) < len(clean_cds) + 200:
            warnings.append(
                WarningMessage(
                    code="CASSETTE_TOO_SHORT",
                    message=(
                        f"Cassette is {len(cassette)} bp but expected at least "
                        f"{len(clean_cds) + 200} bp with flanking elements. "
                        "Promoter or terminator sequences may be incomplete."
                    ),
                )
            )

        return cassette, clean_cds, promoter_seq, feature_order, annotations, cassette_feature_dicts_rebuilt

    # ------------------------------------------------------------------
    # Sequence scanning helpers
    # ------------------------------------------------------------------

    def _find_all_occurrences(self, haystack: str, needle: str) -> List[Tuple[int, int]]:
        """
        Return all non-overlapping exact matches of needle in haystack.
        """
        hits: List[Tuple[int, int]] = []
        if not haystack or not needle:
            return hits
        start = 0
        while True:
            idx = haystack.find(needle, start)
            if idx == -1:
                break
            hits.append((idx, idx + len(needle)))
            start = idx + 1
        return hits

    def _find_backbone_feature_hits_from_insert(
        self,
        backbone: str,
        cassette_feature_dicts: List[Dict],
    ) -> List[BackboneFeatureHit]:
        """
        Scan the backbone directly for the insert features themselves.

        This is critical: even if the GenBank does not annotate CMV/Kozak/polyA/etc.,
        we still find them if the sequence is already present in the backbone.
        """
        hits: List[BackboneFeatureHit] = []

        for feat in cassette_feature_dicts:
            feat_type = self._normalize_feature_type(feat["type"])
            seq = self._clean_seq(feat["sequence"])
            if len(seq) < 6:
                continue

            for start, end in self._find_all_occurrences(backbone, seq):
                hits.append(
                    BackboneFeatureHit(
                        name=feat["name"],
                        type=feat_type,
                        start=start,
                        end=end,
                        sequence=seq,
                        source="sequence_scan",
                    )
                )

        # Deduplicate exact same hit
        uniq: Dict[Tuple[str, str, int, int], BackboneFeatureHit] = {}
        for h in hits:
            uniq[(h.name, h.type, h.start, h.end)] = h
        out = list(uniq.values())

        self._logger.info(
            "[ASSEMBLY] Sequence-scan backbone hits: %s",
            [f"{h.name}:{h.type}[{h.start}:{h.end}]" for h in out],
        )
        return out

    def _merge_backbone_hits(
        self,
        genbank_hits: List[BackboneFeatureHit],
        scan_hits: List[BackboneFeatureHit],
    ) -> List[BackboneFeatureHit]:
        """
        Merge GenBank-derived and sequence-derived feature hits.
        """
        merged: Dict[Tuple[str, str, int, int], BackboneFeatureHit] = {}
        for h in genbank_hits + scan_hits:
            merged[(h.name, h.type, h.start, h.end)] = h
        return sorted(merged.values(), key=lambda x: (x.start, x.end))

    def _choose_replacement_window(
        self,
        backbone: str,
        all_backbone_hits: List[BackboneFeatureHit],
        cassette_feature_dicts: List[Dict],
    ) -> Tuple[int, int, List[BackboneFeatureHit]]:
        """
        Decide which region of the backbone to remove before insertion.

        Strategy:
        - Find all backbone hits that correspond to cassette-like features.
        - Prefer sequence-scan hits because they represent exact insert-content matches.
        - Remove the minimal span covering all duplicated insert features.
        - If nothing is duplicated, fall back to the first exact CDS/promoter match if any.
        - If still nothing, use the MCS-ish fallback backbone_len // 3 insertion site.
        """
        relevant: List[BackboneFeatureHit] = []
        allowed_types = {self._normalize_feature_type(t) for t in _CASSETTE_FEATURE_TYPES}

        for hit in all_backbone_hits:
            if hit.type in allowed_types:
                relevant.append(hit)

        # Prefer exact sequence-scan hits to avoid trusting incomplete GenBank annotation
        scan_hits = [h for h in relevant if h.source == "sequence_scan"]

        if scan_hits:
            cut_start = min(h.start for h in scan_hits)
            cut_end = max(h.end for h in scan_hits)
            return cut_start, cut_end, scan_hits

        if relevant:
            cut_start = min(h.start for h in relevant)
            cut_end = max(h.end for h in relevant)
            return cut_start, cut_end, relevant

        # No duplicated feature found in backbone: insert without replacement
        site = len(backbone) // 3
        return site, site, []

    # ------------------------------------------------------------------
    # Backbone feature rescue / feature-not-found mitigation
    # ------------------------------------------------------------------

    def _infer_backbone_support_warnings(
        self,
        backbone: str,
        warnings: List[WarningMessage],
    ) -> List[WarningMessage]:
        """
        Remove or soften some 'feature not found' situations when the backbone
        obviously already contains them.

        This does NOT magically solve upstream feature resolution, but it avoids
        reporting false absence when the backbone already supplies the function.
        """
        out: List[WarningMessage] = []
        lower_backbone = backbone.upper()

        # Very lightweight heuristics. You can expand later with proper libraries.
        has_sv40_polya_like = "AATAAA" in lower_backbone or "ATTAAA" in lower_backbone
        has_ori_like = "GCGGATAACAATT" in lower_backbone or "CTGTTGACAATTAATCATC" in lower_backbone
        has_kanr_like = "ATGAGGATCGTTTCGCATGATTGAACAAGATGGATTGCACGCA" in lower_backbone

        for w in warnings:
            msg = (w.message or "").lower()

            if "polyadenylation signal" in msg and has_sv40_polya_like:
                out.append(
                    WarningMessage(
                        code="FEATURE_SUPPLIED_BY_BACKBONE",
                        message="Polyadenylation-related sequence appears already present in the backbone.",
                    )
                )
                continue

            if "origin of replication" in msg and has_ori_like:
                out.append(
                    WarningMessage(
                        code="FEATURE_SUPPLIED_BY_BACKBONE",
                        message="Origin-of-replication-related sequence appears already present in the backbone.",
                    )
                )
                continue

            if "selection markers" in msg and has_kanr_like:
                out.append(
                    WarningMessage(
                        code="FEATURE_SUPPLIED_BY_BACKBONE",
                        message="Selection-marker-related sequence appears already present in the backbone.",
                    )
                )
                continue

            out.append(w)

        return out

    # ------------------------------------------------------------------
    # Duplication check
    # ------------------------------------------------------------------

    def _duplication_warnings(
        self,
        assembled: str,
        cassette_feature_dicts: List[Dict],
    ) -> List[WarningMessage]:
        """
        Check if any feature appears multiple times in the assembled plasmid.
        """
        warnings: List[WarningMessage] = []

        for feat in cassette_feature_dicts:
            seq = self._clean_seq(feat["sequence"])
            if len(seq) < 20:
                continue

            count = assembled.count(seq)
            if count > 1:
                warnings.append(
                    WarningMessage(
                        code="FEATURE_DUPLICATED_IN_ASSEMBLY",
                        message=(
                            f"Feature '{feat['name']}' appears {count}x in the assembled plasmid. "
                            "Backbone replacement may have been incomplete."
                        ),
                    )
                )

        return warnings

    # ------------------------------------------------------------------
    # Assembly methods
    # ------------------------------------------------------------------

    def _assemble_gibson(
        self, backbone: str, cassette: str, cut_start: int, cut_end: int
    ) -> Tuple[str, List[AssemblyFragment], List[str]]:
        """
        Gibson assembly: design overlapping sequences for isothermal assembly.

        IT (user-facing explanation):
        Nel caso della Gibson Assembly, progetta il frammento aggiungendo alle estremità
        sinistra e destra regioni di omologia (~20–40 bp) che corrispondono alle sequenze
        adiacenti sul backbone. Poi amplifica il frammento tramite PCR usando primer che
        incorporano queste code di omologia, così da ottenere sovrapposizioni accurate e
        facilitare l’assemblaggio one‑pot.
        """
        hom_len = 30

        left_overlap = backbone[max(0, cut_start - hom_len):cut_start]
        right_overlap = backbone[cut_end:cut_end + hom_len]

        assembled = backbone[:cut_start] + cassette + backbone[cut_end:]

        fwd_primer = f"Fwd: 5'-{left_overlap[-20:]}{cassette[:20]}-3'"
        rev_primer = f"Rev: 5'-{_reverse_complement(right_overlap[:20])}{_reverse_complement(cassette[-20:])}-3'"

        fragments = [
            AssemblyFragment(
                name="backbone_left_arm",
                sequence=backbone[:cut_start],
                left_homology=None,
                right_homology=left_overlap or None,
            ),
            AssemblyFragment(
                name="expression_cassette",
                sequence=cassette,
                left_homology=left_overlap or None,
                right_homology=right_overlap or None,
            ),
            AssemblyFragment(
                name="backbone_right_arm",
                sequence=backbone[cut_end:],
                left_homology=right_overlap or None,
                right_homology=None,
            ),
        ]
        return assembled, fragments, [fwd_primer, rev_primer]

    def _assemble_golden_gate(
        self, backbone: str, cassette: str, cut_start: int, cut_end: int
    ) -> Tuple[str, List[AssemblyFragment], List[str]]:
        """
        Golden Gate assembly: Type IIS restriction site-based modular cloning.
        """
        assembled = backbone[:cut_start] + cassette + backbone[cut_end:]
        fragments = [
            AssemblyFragment(name="backbone_fragment", sequence=backbone, left_homology=None, right_homology=None),
            AssemblyFragment(name="expression_cassette", sequence=cassette, left_homology=None, right_homology=None),
        ]
        primers = [
            "Design primers with flanking BsaI or BsmBI type IIS restriction sites.",
            "Ensure 4 bp overhangs are unique and directional for ordered ligation.",
            "Verify absence of internal BsaI/BsmBI sites in insert and backbone.",
        ]
        return assembled, fragments, primers

    def _assemble_restriction(
        self, backbone: str, cassette: str, cut_start: int, cut_end: int
    ) -> Tuple[str, List[AssemblyFragment], List[str]]:
        """
        Restriction enzyme digestion and ligation.
        """
        assembled = backbone[:cut_start] + cassette + backbone[cut_end:]
        fragments = [
            AssemblyFragment(name="backbone_fragment", sequence=backbone, left_homology=None, right_homology=None),
            AssemblyFragment(name="expression_cassette", sequence=cassette, left_homology=None, right_homology=None),
        ]
        primers = [
            "Add compatible restriction enzyme sites to both insert ends.",
            "Digest backbone and insert with matching enzymes (e.g. EcoRI/XhoI).",
            "Verify that no internal forbidden restriction sites are present in the insert.",
        ]
        return assembled, fragments, primers

    def _assemble_synthesis(
        self, backbone: str, cassette: str, cut_start: int, cut_end: int
    ) -> Tuple[str, List[AssemblyFragment], List[str]]:
        """
        Chemical synthesis: submit sequences for de novo synthesis.
        """
        assembled = backbone[:cut_start] + cassette + backbone[cut_end:]
        fragments = [
            AssemblyFragment(name="full_backbone_fragment", sequence=backbone, left_homology=None, right_homology=None),
            AssemblyFragment(name="synthetic_insert_fragment", sequence=cassette, left_homology=None, right_homology=None),
        ]
        primers = ["No PCR primers required; submit for chemical synthesis."]
        return assembled, fragments, primers

    # ------------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------------

    def run(self, inp: AssemblyInput) -> AssemblyOutput:
        """
        Main entry point for assembly planning.

        Flow:
        1. Build initial cassette with all features
        2. Parse backbone for existing features (GenBank + sequence scan)
        3. Filter cassette to remove duplicate features already in backbone
        4. Rebuild cassette with only non-duplicate features
        5. Choose replacement window in backbone
        6. Apply assembly method
        7. Final validation and warning refinement
        """
        self._logger.info("[ASSEMBLY] INPUT backbone=%s method=%s", inp.backbone_name, inp.assembly_method)

        warnings: List[WarningMessage] = []

        # Step 1: Build initial cassette
        cassette, cds_seq, promoter_seq, feature_order, cassette_annotations, cassette_feature_dicts = self._build_cassette(
            inp.features, inp.cds_sequence, inp.gene_symbol, warnings
        )

        backbone = self._clean_seq(inp.backbone_sequence)
        method = self._pick_method(inp.assembly_method)

        if not backbone:
            warnings.append(WarningMessage(code="EMPTY_BACKBONE", message="Backbone sequence is empty."))
        if not cassette:
            warnings.append(WarningMessage(code="EMPTY_INSERT", message="Construct sequence is empty."))

        # Step 2: Parse existing GenBank annotations
        genbank_hits = self._parse_genbank_features(inp.backbone_genbank)

        # Step 2b: Sequence-level scan for insert features already present in backbone
        scan_hits = self._find_backbone_feature_hits_from_insert(backbone, cassette_feature_dicts)

        # Step 2c: Merge both sources
        all_backbone_hits = self._merge_backbone_hits(genbank_hits, scan_hits)

        self._logger.info(
            "[ASSEMBLY] Backbone feature universe: %s",
            [f"{h.name}:{h.type}[{h.start}:{h.end}]/{h.source}" for h in all_backbone_hits],
        )

        # Step 3: Filter out cassette features already in backbone
        filtered_cassette_feature_dicts, filter_warnings = self._filter_cassette_features(
            cassette_feature_dicts,
            backbone,
            all_backbone_hits,
        )
        warnings.extend(filter_warnings)

        # Step 4: Rebuild cassette with only the filtered features
        cassette, cds_seq, promoter_seq, feature_order, cassette_annotations, cassette_feature_dicts = self._rebuild_cassette_from_filtered(
            inp.features,
            filtered_cassette_feature_dicts,
            inp.cds_sequence,
            inp.gene_symbol,
            warnings,
        )

        # Step 5: Choose replacement window
        cut_start, cut_end, replaced_hits = self._choose_replacement_window(
            backbone=backbone,
            all_backbone_hits=all_backbone_hits,
            cassette_feature_dicts=cassette_feature_dicts,
        )

        if replaced_hits:
            replaced_names = [f"{h.name}:{h.type}" for h in replaced_hits]
            self._logger.info(
                "[ASSEMBLY] Replacing existing backbone region [%d:%d] containing: %s",
                cut_start,
                cut_end,
                replaced_names,
            )
            warnings.append(
                WarningMessage(
                    code="BACKBONE_FEATURES_REPLACED",
                    message=(
                        f"Existing backbone feature region [{cut_start}:{cut_end}] was removed and replaced with the "
                        f"new cassette. Replaced features: {', '.join(replaced_names)}."
                    ),
                )
            )
        else:
            self._logger.info(
                "[ASSEMBLY] No existing insert-equivalent feature detected in backbone; "
                "using fallback insertion site at %d.",
                cut_start,
            )

        # Step 6: Dispatch to assembly method
        dispatch = {
            "Gibson": self._assemble_gibson,
            "Golden Gate": self._assemble_golden_gate,
            "Restriction": self._assemble_restriction,
        }
        assembler = dispatch.get(method, self._assemble_synthesis)
        assembled, fragments, primers = assembler(backbone, cassette, cut_start, cut_end)

        # Step 7: Re-check unresolved-feature warnings against actual backbone content
        warnings = self._infer_backbone_support_warnings(backbone, warnings)

        # Step 7b: Post-assembly duplication check
        warnings.extend(self._duplication_warnings(assembled, cassette_feature_dicts))

        junction_offset = cut_start + 1

        out = AssemblyOutput(
            assembly_method=method,
            fragments=fragments,
            primer_requirements=primers,
            assembled_sequence=assembled,
            junction_offset=junction_offset,
            cassette_sequence=cassette,
            feature_order=feature_order,
            cassette_annotations=cassette_annotations,
            warnings=warnings,
        )

        self._logger.info("[ASSEMBLY] OUTPUT sequence_len=%d warnings=%d", len(assembled), len(warnings))
        return out