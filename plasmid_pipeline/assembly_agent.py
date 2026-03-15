from __future__ import annotations

import io
import logging
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from logging_utils import get_conversation_logger

from .models import AssemblyInput, AssemblyOutput, AssemblyFragment, Annotation, WarningMessage

_RC_TABLE = str.maketrans("ATGCN", "TACGN")

logger = logging.getLogger(__name__)


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


# Feature types considered "expression cassette" features — if they exist in
# both backbone and cassette they must be deduplicated.
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
}


class AssemblyAgent:
    """
    Biologically-aware, GenBank-driven assembly planner.

    For each assembly:
    1.  Parse the backbone GenBank record (if present) to extract annotated features.
    2.  Compare backbone features against the expression cassette features
        (from ConstructOutput.annotations) to find conflicts.
    3.  Derive cut sites automatically: left arm ends just before the first
        conflicting backbone feature; right arm starts just after the last.
    4.  Fall back to backbone_len // 3 when no GenBank record is available or
        no conflicts are found.
    5.  Validate the final sequence — no feature sequence should appear twice.
    """

    def __init__(self) -> None:
        self._logger = get_conversation_logger()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _clean_seq(self, seq: str) -> str:
        return "".join(ch for ch in (seq or "").upper() if ch in {"A", "T", "G", "C", "N"})

    def _classify_host(self, expression_host: str) -> str:
        host = (expression_host or "").lower()
        for kw in ("hek", "cho", "cos", "vero", "hela", "293", "mammalian", "human", "mouse", "rat"):
            if kw in host:
                return "mammalian"
        for kw in ("e.coli", "ecoli", "bl21", "dh5", "bacteria", "bacterial", "prokaryot"):
            if kw in host:
                return "bacterial"
        return "unknown"

    def _pick_method(self, pref: str) -> str:
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

    # ------------------------------------------------------------------
    # GenBank parsing
    # ------------------------------------------------------------------

    def _parse_genbank_features(self, genbank_text: str) -> List[Dict]:
        """
        Parse raw GenBank text and return a list of dicts:
            {name, type, start, end, sequence}
        for every annotated feature (excluding 'source' and 'gene' loci).
        """
        features: List[Dict] = []
        try:
            record: SeqRecord = SeqIO.read(io.StringIO(genbank_text), "genbank")
        except Exception as exc:
            self._logger.warning("[ASSEMBLY] GenBank parse failed: %s", exc)
            return features

        full_seq = str(record.seq).upper()

        for feat in record.features:
            if feat.type in ("source", "gene"):
                continue

            name = (
                feat.qualifiers.get("label", [""])[0]
                or feat.qualifiers.get("gene", [""])[0]
                or feat.qualifiers.get("product", [""])[0]
                or feat.type
            )

            try:
                start = int(feat.location.start)
                end   = int(feat.location.end)
                seq   = full_seq[start:end]
            except Exception:
                continue

            features.append({
                "name":  name,
                "type":  feat.type,
                "start": start,
                "end":   end,
                "sequence": seq,
            })

        self._logger.info(
            "[ASSEMBLY] Backbone GenBank features: %s",
            [f["name"] for f in features],
        )
        return features

    # ------------------------------------------------------------------
    # Conflict detection & cut-site derivation
    # ------------------------------------------------------------------

    def _build_cassette(
        self,
        features,
        cds_sequence: str,
        gene_symbol: str,
        warnings: List[WarningMessage],
    ) -> Tuple[str, str, str, List[str], List[Annotation]]:
        """
        Build the expression cassette:
            promoter + kozak + [n_tags] + cds + [unspecified_tags] + [c_tags] + terminator

        Returns (cassette, clean_cds, promoter_seq, feature_order, annotations).
        CDS validation warnings are appended to the shared warnings list.
        """
        promoters   = [f for f in features if f.type == "promoter"]
        kozak_feats = [f for f in features if f.type == "kozak"]
        tag_feats   = [f for f in features if f.type == "tag"]
        terminators = [f for f in features if f.type == "terminator"]

        n_tags          = [t for t in tag_feats if (t.position_hint or "").upper() == "N"]
        c_tags          = [t for t in tag_feats if (t.position_hint or "").upper() == "C"]
        unspecified_tags = [t for t in tag_feats if (t.position_hint or "").upper() not in ("N", "C")]

        if unspecified_tags:
            warnings.append(WarningMessage(
                code="TAG_POSITION_UNSPECIFIED",
                message=(
                    "One or more tag features have no position_hint. "
                    "Unspecified tags will be placed after the CDS by default."
                ),
            ))

        parts: List[str] = []
        feature_order: List[str] = []
        annotations: List[Annotation] = []
        cursor = 1

        def _append(name: str, seq: str, feat_type: str, label: str) -> None:
            nonlocal cursor
            clean = self._clean_seq(seq)
            if not clean:
                return
            parts.append(clean)
            annotations.append(Annotation(name=name, start=cursor, end=cursor + len(clean) - 1, type=feat_type))
            feature_order.append(label)
            cursor += len(clean)

        # Promoter
        if promoters:
            p = promoters[0]
            _append(p.name, p.sequence, "promoter", f"{p.name}_promoter")
            self._logger.info("[ASSEMBLY] promoter (%s): %d bp", p.name, len(self._clean_seq(p.sequence)))
        else:
            warnings.append(WarningMessage(
                code="NO_PROMOTER",
                message="No promoter feature resolved; construct may not be expressed.",
            ))
            self._logger.info("[ASSEMBLY] promoter: 0 bp (not resolved)")

        # Kozak
        if kozak_feats:
            k = kozak_feats[0]
            _append("Kozak", k.sequence, "misc_feature", "Kozak")
            self._logger.info("[ASSEMBLY] kozak: %d bp", len(self._clean_seq(k.sequence)))
        else:
            self._logger.info("[ASSEMBLY] kozak: 0 bp")

        # N-terminal tags
        for tag in n_tags:
            _append(f"{tag.name}_tag", tag.sequence, "misc_feature", f"{tag.name}_N_tag")

        # CDS — validate then append
        clean_cds = self._clean_seq(cds_sequence)
        if not clean_cds:
            warnings.append(WarningMessage(code="NO_CDS", message="No CDS sequence available; construct will be incomplete."))
        else:
            if len(clean_cds) % 3 != 0:
                warnings.append(WarningMessage(code="CDS_NOT_MULTIPLE_OF_3", message="CDS length is not divisible by 3."))
            if not clean_cds.startswith("ATG"):
                warnings.append(WarningMessage(code="CDS_NO_START_CODON", message="CDS does not start with ATG."))
            if clean_cds[-3:] not in {"TAA", "TAG", "TGA"}:
                warnings.append(WarningMessage(code="CDS_NO_STOP_CODON", message="CDS does not end with a standard stop codon."))
            _append(f"{gene_symbol}_CDS", clean_cds, "CDS", f"{gene_symbol}_CDS")
            self._logger.info("[ASSEMBLY] cds: %d bp", len(clean_cds))

        # Unspecified tags (after CDS)
        for tag in unspecified_tags:
            _append(f"{tag.name}_tag", tag.sequence, "misc_feature", f"{tag.name}_tag")

        # C-terminal tags
        for tag in c_tags:
            _append(f"{tag.name}_tag", tag.sequence, "misc_feature", f"{tag.name}_C_tag")

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
            warnings.append(WarningMessage(
                code="CASSETTE_TOO_SHORT",
                message=(
                    f"Cassette is {len(cassette)} bp but expected at least "
                    f"{len(clean_cds) + 200} bp with flanking elements. "
                    "Promoter or terminator sequences may be incomplete."
                ),
            ))

        return cassette, clean_cds, promoter_seq, feature_order, annotations

    def _cassette_feature_seqs(self, inp: AssemblyInput) -> List[str]:
        """
        Return a list of cleaned sequences for every feature in the cassette.
        """
        seqs: List[str] = []
        for f in inp.features:
            seq = self._clean_seq(f.sequence)
            if len(seq) >= 20:
                seqs.append(seq)
        return seqs

    def _find_cut_sites(
        self,
        backbone: str,
        backbone_features: List[Dict],
        cassette_feature_seqs: List[str],
    ) -> Tuple[int, List[str], List[str]]:
        """
        Returns (insertion_site, conflicting_feature_names, all_backbone_feature_names).

        Strategy:
        - A backbone feature "conflicts" if its sequence appears in the cassette.
        - Cut site = span covering all conflicting features:
            left  arm ends at  min(conflict_starts)
            right arm starts at max(conflict_ends)
        - Falls back to backbone_len // 3 when no conflicts found.
        """
        all_names = [f["name"] for f in backbone_features]
        conflicting: List[Dict] = []

        for bf in backbone_features:
            if bf["type"] not in _CASSETTE_FEATURE_TYPES:
                continue
            bf_seq = bf["sequence"]
            if len(bf_seq) < 20:
                continue
            for cf_seq in cassette_feature_seqs:
                # Check a 40-bp probe from each end to allow partial matches
                probe = bf_seq[:40]
                if probe and probe in cf_seq:
                    conflicting.append(bf)
                    break
                probe2 = bf_seq[-40:]
                if probe2 and probe2 in cf_seq:
                    conflicting.append(bf)
                    break

        conflict_names = [f["name"] for f in conflicting]

        if conflicting:
            cut_start = min(f["start"] for f in conflicting)
            cut_end   = max(f["end"]   for f in conflicting)
            self._logger.info(
                "[ASSEMBLY] Conflicting backbone features: %s — cutting [%d:%d]",
                conflict_names, cut_start, cut_end,
            )
            # insertion_site = start of the removed span; right arm picks up at cut_end
            return cut_start, conflict_names, all_names
        else:
            site = len(backbone) // 3
            self._logger.info(
                "[ASSEMBLY] No conflicting features found; using MCS fallback site=%d", site
            )
            return site, [], all_names

    # ------------------------------------------------------------------
    # Post-assembly duplication check
    # ------------------------------------------------------------------

    def _duplication_warnings(
        self,
        assembled: str,
        cassette_feature_seqs: List[str],
    ) -> List[WarningMessage]:
        warnings: List[WarningMessage] = []
        for seq in cassette_feature_seqs:
            if len(seq) < 40:
                continue
            probe = seq[:40]
            count = assembled.count(probe)
            if count > 1:
                warnings.append(
                    WarningMessage(
                        code="FEATURE_DUPLICATED_IN_ASSEMBLY",
                        message=(
                            f"A cassette feature sequence (probe: {probe[:20]}…) "
                            f"appears {count}x in the assembled plasmid. "
                            "Check that the backbone was properly linearised."
                        ),
                    )
                )
        return warnings

    # ------------------------------------------------------------------
    # Assembly methods
    # ------------------------------------------------------------------

    def _assemble_gibson(
        self, backbone: str, cassette: str, site: int, cut_end: Optional[int] = None
    ) -> Tuple[str, List[AssemblyFragment], List[str]]:
        hom_len = 30
        left_arm  = backbone[max(0, site - hom_len) : site]
        end       = cut_end if cut_end is not None else site
        right_arm = backbone[end : end + hom_len]

        assembled = backbone[:site] + cassette + backbone[end:]

        fwd_primer = f"Fwd: 5'-{left_arm[-20:]}{cassette[:20]}-3'"
        rev_primer = f"Rev: 5'-{_reverse_complement(right_arm[:20])}{_reverse_complement(cassette[-20:])}-3'"

        fragments: List[AssemblyFragment] = [
            AssemblyFragment(
                name="backbone_left_arm",
                sequence=backbone[:site],
                left_homology=None,
                right_homology=left_arm or None,
            ),
            AssemblyFragment(
                name="expression_cassette",
                sequence=cassette,
                left_homology=left_arm or None,
                right_homology=right_arm or None,
            ),
            AssemblyFragment(
                name="backbone_right_arm",
                sequence=backbone[end:],
                left_homology=right_arm or None,
                right_homology=None,
            ),
        ]
        return assembled, fragments, [fwd_primer, rev_primer]

    def _assemble_golden_gate(
        self, backbone: str, cassette: str, site: int, cut_end: Optional[int] = None
    ) -> Tuple[str, List[AssemblyFragment], List[str]]:
        end = cut_end if cut_end is not None else site
        assembled = backbone[:site] + cassette + backbone[end:]
        fragments: List[AssemblyFragment] = [
            AssemblyFragment(name="backbone_fragment",    sequence=backbone,  left_homology=None, right_homology=None),
            AssemblyFragment(name="expression_cassette",  sequence=cassette,  left_homology=None, right_homology=None),
        ]
        primers = [
            "Design primers with flanking BsaI or BsmBI type IIS restriction sites.",
            "Ensure 4 bp overhangs are unique and directional for ordered ligation.",
            "Verify absence of internal BsaI/BsmBI sites in insert and backbone.",
        ]
        return assembled, fragments, primers

    def _assemble_restriction(
        self, backbone: str, cassette: str, site: int, cut_end: Optional[int] = None
    ) -> Tuple[str, List[AssemblyFragment], List[str]]:
        end = cut_end if cut_end is not None else site
        assembled = backbone[:site] + cassette + backbone[end:]
        fragments: List[AssemblyFragment] = [
            AssemblyFragment(name="backbone_fragment",    sequence=backbone,  left_homology=None, right_homology=None),
            AssemblyFragment(name="expression_cassette",  sequence=cassette,  left_homology=None, right_homology=None),
        ]
        primers = [
            "Add compatible restriction enzyme sites to both insert ends.",
            "Digest backbone and insert with matching enzymes (e.g. EcoRI/XhoI).",
            "Verify that no internal forbidden restriction sites are present in the insert.",
        ]
        return assembled, fragments, primers

    def _assemble_synthesis(
        self, backbone: str, cassette: str, site: int, cut_end: Optional[int] = None
    ) -> Tuple[str, List[AssemblyFragment], List[str]]:
        end = cut_end if cut_end is not None else site
        assembled = backbone[:site] + cassette + backbone[end:]
        fragments: List[AssemblyFragment] = [
            AssemblyFragment(name="full_backbone_fragment",   sequence=backbone, left_homology=None, right_homology=None),
            AssemblyFragment(name="synthetic_insert_fragment", sequence=cassette, left_homology=None, right_homology=None),
        ]
        primers = ["No PCR primers required; submit for chemical synthesis."]
        return assembled, fragments, primers

    # ------------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------------

    def run(self, inp: AssemblyInput) -> AssemblyOutput:
        self._logger.info("[ASSEMBLY] INPUT backbone=%s method=%s", inp.backbone_name, inp.assembly_method)

        warnings: List[WarningMessage] = []
        cassette, cds_seq, promoter_seq, feature_order, cassette_annotations = self._build_cassette(
            inp.features, inp.cds_sequence, inp.gene_symbol, warnings
        )

        backbone = self._clean_seq(inp.backbone_sequence)
        method   = self._pick_method(inp.assembly_method)

        if not backbone:
            warnings.append(WarningMessage(code="EMPTY_BACKBONE", message="Backbone sequence is empty."))
        if not cassette:
            warnings.append(WarningMessage(code="EMPTY_INSERT",   message="Construct sequence is empty."))

        # --- Parse backbone GenBank and find cut sites ---
        backbone_features: List[Dict] = []
        if inp.backbone_genbank:
            backbone_features = self._parse_genbank_features(inp.backbone_genbank)

        cassette_feat_seqs = self._cassette_feature_seqs(inp)

        cut_start, conflict_names, all_bb_names = self._find_cut_sites(
            backbone, backbone_features, cassette_feat_seqs
        )

        # cut_end is only meaningful when there are conflicts
        cut_end: Optional[int] = None
        if conflict_names and backbone_features:
            conflicting = [f for f in backbone_features if f["name"] in conflict_names]
            cut_end = max(f["end"] for f in conflicting)

        # Log what was found / removed
        self._logger.info("[ASSEMBLY] All backbone features: %s", all_bb_names)
        if conflict_names:
            self._logger.info(
                "[ASSEMBLY] Removing conflicting backbone features to make room for cassette: %s",
                conflict_names,
            )
            warnings.append(
                WarningMessage(
                    code="BACKBONE_FEATURES_REMOVED",
                    message=(
                        f"The following backbone features overlap with the expression cassette "
                        f"and were removed during linearisation: {', '.join(conflict_names)}."
                    ),
                )
            )
        else:
            self._logger.info(
                "[ASSEMBLY] No backbone features conflict with the cassette; using fallback insertion site."
            )

        # --- Dispatch to assembly method ---
        dispatch = {
            "Gibson":       self._assemble_gibson,
            "Golden Gate":  self._assemble_golden_gate,
            "Restriction":  self._assemble_restriction,
        }
        assembler = dispatch.get(method, self._assemble_synthesis)
        assembled, fragments, primers = assembler(backbone, cassette, cut_start, cut_end)

        # --- Post-assembly duplication check ---
        warnings.extend(self._duplication_warnings(assembled, cassette_feat_seqs))

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
