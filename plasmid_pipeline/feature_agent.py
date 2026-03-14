from __future__ import annotations

from typing import List, Optional

from logging_utils import get_conversation_logger

from .models import FeatureInput, FeatureOutput, ResolvedFeature, WarningMessage


class FeatureAgent:
    """
    Resolves promoter / kozak / terminator features with strict biological checks.

    Key fixes:
    - Rejects ultra-short fake promoters.
    - Separates promoter from Kozak.
    - Adds a real default mammalian polyA if terminator/polyA is missing.
    """

    # Minimal practical promoter lengths to avoid misclassifying Kozak/start context.
    MIN_PROMOTER_LENGTH = 80

    # Canonical default parts.
    CMV_PROMOTER_SEQ = (
        "CGCAAATGGGCGGTAGGCGTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGG"
        "CGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGT"
        "GCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAGCG"
    )

    KOZAK_SEQ = "GCCACC"

    BGH_POLYA_SEQ = (
        "AATAAA"
        "GATCTTTATTTTCATTAGATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGCGATAAGGATCC"
    )

    SV40_POLYA_SEQ = (
        "AATAAA"
        "GCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTATCATG"
    )

    def __init__(self) -> None:
        self._logger = get_conversation_logger()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _clean_sequence(self, seq: Optional[str]) -> str:
        if not seq:
            return ""
        seq = seq.upper().replace("\n", "").replace("\r", "").replace(" ", "")
        return "".join(ch for ch in seq if ch in {"A", "T", "G", "C", "N"})

    def _looks_like_fake_promoter(self, seq: str) -> bool:
        """
        Reject short sequences that are really Kozak/start-context fragments.
        """
        seq = self._clean_sequence(seq)
        if len(seq) < self.MIN_PROMOTER_LENGTH:
            return True

        # Strong sign of start-context masquerading as promoter
        if seq in {"GCCACC", "GCCGCCACCATGG", "CCACCATG", "ACCATG"}:
            return True

        return False

    def _default_promoter_for_input(self, inp: FeatureInput) -> Optional[ResolvedFeature]:
        promoter_name = (inp.promoter or "").upper()
        host = (inp.expression_host or "").lower()

        if promoter_name == "CMV" and any(x in host for x in ["hek", "293", "cho", "mamm"]):
            return ResolvedFeature(
                name="CMV",
                type="promoter",
                sequence=self.CMV_PROMOTER_SEQ,
                position_hint=None,
            )

        return None

    def _default_kozak_for_input(self, inp: FeatureInput) -> Optional[ResolvedFeature]:
        host = (inp.expression_host or "").lower()
        if any(x in host for x in ["hek", "293", "cho", "mamm"]):
            return ResolvedFeature(
                name="Kozak",
                type="kozak",
                sequence=self.KOZAK_SEQ,
                position_hint=None,
            )
        return None

    def _default_terminator_for_input(self, inp: FeatureInput) -> Optional[ResolvedFeature]:
        host = (inp.expression_host or "").lower()
        if any(x in host for x in ["hek", "293", "cho", "mamm"]):
            return ResolvedFeature(
                name="bGH_polyA",
                type="terminator",
                sequence=self.BGH_POLYA_SEQ,
                position_hint=None,
            )
        return None

    def _normalize_feature(self, feat: ResolvedFeature, warnings: List[WarningMessage]) -> Optional[ResolvedFeature]:
        seq = self._clean_sequence(feat.sequence)

        if feat.type == "promoter":
            if self._looks_like_fake_promoter(seq):
                warnings.append(
                    WarningMessage(
                        code="REJECTED_FAKE_OR_TOO_SHORT_PROMOTER",
                        message=(
                            f"Rejected promoter candidate '{feat.name}' because its sequence is too short "
                            f"({len(seq)} bp) or looks like a Kozak/start-context fragment."
                        ),
                    )
                )
                return None

        if feat.type == "kozak":
            if seq != self.KOZAK_SEQ:
                warnings.append(
                    WarningMessage(
                        code="NONCANONICAL_KOZAK",
                        message=f"Kozak sequence '{feat.name}' is non-canonical; keeping as provided.",
                    )
                )

        return ResolvedFeature(
            name=feat.name,
            type=feat.type,
            sequence=seq,
            position_hint=feat.position_hint,
        )

    # ------------------------------------------------------------------
    # Main
    # ------------------------------------------------------------------

    def run(self, inp: FeatureInput) -> FeatureOutput:
        self._logger.info("[FEATURE] INPUT %s", inp.model_dump())
        warnings: List[WarningMessage] = []
        out_features: List[ResolvedFeature] = []

        # 1) Normalize user/upstream-provided features first.
        incoming = list(inp.features or [])
        for feat in incoming:
            normalized = self._normalize_feature(feat, warnings)
            if normalized is not None:
                out_features.append(normalized)

        existing_types = {f.type for f in out_features}
        existing_names_upper = {f.name.upper() for f in out_features}

        # 2) Promoter
        if "promoter" not in existing_types and (inp.promoter or "").strip():
            promoter = self._default_promoter_for_input(inp)
            if promoter:
                out_features.append(promoter)
            else:
                warnings.append(
                    WarningMessage(
                        code="PROMOTER_NOT_RESOLVED",
                        message=f"Could not resolve a valid promoter sequence for requested promoter '{inp.promoter}'.",
                    )
                )

        # 3) Kozak
        if "kozak" not in existing_types:
            kozak = self._default_kozak_for_input(inp)
            if kozak:
                out_features.append(kozak)

        # 4) Terminator / polyA
        has_terminator = any(f.type == "terminator" for f in out_features)
        if not has_terminator:
            default_term = self._default_terminator_for_input(inp)
            if default_term:
                out_features.append(default_term)
                warnings.append(
                    WarningMessage(
                        code="DEFAULT_TERMINATOR_ADDED",
                        message=(
                            f"No terminator/polyA was provided. Added default mammalian terminator '{default_term.name}'."
                        ),
                    )
                )

        # 5) Sanity checks
        promoter_feats = [f for f in out_features if f.type == "promoter"]
        if len(promoter_feats) > 1:
            warnings.append(
                WarningMessage(
                    code="MULTIPLE_PROMOTERS_RESOLVED",
                    message="Multiple promoters were resolved. Downstream construct assembly should choose only one.",
                )
            )

        kozak_feats = [f for f in out_features if f.type == "kozak"]
        if len(kozak_feats) > 1:
            warnings.append(
                WarningMessage(
                    code="MULTIPLE_KOZAKS_RESOLVED",
                    message="Multiple Kozak sequences were resolved. Downstream construct assembly should choose only one.",
                )
            )

        terminator_feats = [f for f in out_features if f.type == "terminator"]
        if len(terminator_feats) > 1:
            warnings.append(
                WarningMessage(
                    code="MULTIPLE_TERMINATORS_RESOLVED",
                    message="Multiple terminator/polyA features were resolved. Downstream construct assembly should choose only one.",
                )
            )

        out = FeatureOutput(features=out_features, warnings=warnings)
        self._logger.info(
            "[FEATURE] OUTPUT features=%s warnings=%s",
            [f.model_dump() for f in out.features],
            [w.model_dump() for w in out.warnings],
        )
        return out