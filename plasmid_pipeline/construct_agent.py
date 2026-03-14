from __future__ import annotations

from typing import List

from logging_utils import get_conversation_logger

from .models import (
    ConstructInput,
    ConstructOutput,
    Annotation,
    ResolvedFeature,
    WarningMessage,
)


class ConstructAgent:
    """
    Builds a linear insert cassette, now backbone-aware.

    Rules:
    - If the backbone already appears to contain the requested promoter, do not add it again.
    - If the backbone looks like a loaded vector, default to a reduced insert strategy.
    - Keep tag/CDS ordering biologically sensible.
    """

    def __init__(self) -> None:
        self._logger = get_conversation_logger()

    def _find_feature(self, feats: List[ResolvedFeature], feat_type: str) -> List[ResolvedFeature]:
        return [f for f in feats if f.type == feat_type]

    def _clean_seq(self, seq: str) -> str:
        return "".join(ch for ch in (seq or "").upper() if ch in {"A", "T", "G", "C", "N"})

    def _cds_warnings(self, cds: str) -> List[WarningMessage]:
        warnings: List[WarningMessage] = []
        if not cds:
            warnings.append(
                WarningMessage(
                    code="NO_CDS",
                    message="No CDS sequence available; construct will be incomplete.",
                )
            )
            return warnings

        if len(cds) % 3 != 0:
            warnings.append(
                WarningMessage(
                    code="CDS_NOT_MULTIPLE_OF_3",
                    message="CDS length is not divisible by 3.",
                )
            )

        if not cds.startswith("ATG"):
            warnings.append(
                WarningMessage(
                    code="CDS_NO_START_CODON",
                    message="CDS does not start with ATG.",
                )
            )

        stop_codons = {"TAA", "TAG", "TGA"}
        if len(cds) >= 3 and cds[-3:] not in stop_codons:
            warnings.append(
                WarningMessage(
                    code="CDS_NO_STOP_CODON",
                    message="CDS does not end with a standard stop codon.",
                )
            )

        return warnings

    def _split_tags(self, feats: List[ResolvedFeature]) -> tuple[List[ResolvedFeature], List[ResolvedFeature], List[ResolvedFeature]]:
        n_tags: List[ResolvedFeature] = []
        c_tags: List[ResolvedFeature] = []
        unspecified: List[ResolvedFeature] = []

        for tag in feats:
            hint = (tag.position_hint or "").strip().upper()
            if hint == "N":
                n_tags.append(tag)
            elif hint == "C":
                c_tags.append(tag)
            else:
                unspecified.append(tag)

        return n_tags, c_tags, unspecified

    def _append_feature(
        self,
        *,
        sequence_parts: List[str],
        annotations: List[Annotation],
        feature_order: List[str],
        cursor: int,
        name: str,
        seq: str,
        feat_type: str,
        order_label: str | None = None,
    ) -> int:
        clean = self._clean_seq(seq)
        if not clean:
            return cursor

        sequence_parts.append(clean)
        start = cursor
        end = cursor + len(clean) - 1

        annotations.append(
            Annotation(
                name=name,
                start=start,
                end=end,
                type=feat_type,
            )
        )
        feature_order.append(order_label or name)
        return end + 1

    def run(self, inp: ConstructInput) -> ConstructOutput:
        self._logger.info("[CONSTRUCT] INPUT %s", inp.model_dump())

        sequence_parts: List[str] = []
        annotations: List[Annotation] = []
        feature_order: List[str] = []
        warnings: List[WarningMessage] = []

        cursor = 1
        host_system = (inp.expression.host_system or "").lower()

        promoters = self._find_feature(inp.features.features, "promoter")
        kozak_feats = self._find_feature(inp.features.features, "kozak")
        tag_feats = self._find_feature(inp.features.features, "tag")
        terminators = self._find_feature(inp.features.features, "terminator")

        n_tags, c_tags, unspecified_tags = self._split_tags(tag_feats)

        requested_promoter_names = [p.name.upper() for p in promoters]
        backbone_promoters = [p.upper() for p in inp.backbone.backbone_promoters]

        include_promoter = True
        include_terminator = True
        insertion_mode: str = "full_cassette"

        if inp.backbone.is_loaded_vector:
            warnings.append(
                WarningMessage(
                    code="LOADED_BACKBONE_CONTEXT",
                    message=(
                        "Selected backbone looks like a loaded vector. "
                        "ConstructAgent will avoid adding a second full cassette when possible."
                    ),
                )
            )
            insertion_mode = "cds_plus_tags"

        # Heuristic: if backbone name clearly indicates an expression vector (e.g. pCMV)
        # and a promoter is requested, avoid duplicating the promoter cassette.
        name_lower = inp.backbone.backbone_name.lower()
        if any(x in name_lower for x in ["expression vector", "pcmv", "ef1", "cag", "pgk"]):
            if promoters:
                include_promoter = False
                insertion_mode = "cds_plus_tags"
                warnings.append(
                    WarningMessage(
                        code="BACKBONE_EXPRESSION_VECTOR_HEURISTIC",
                        message=(
                            "Backbone name suggests it is an expression vector; "
                            "requested promoter will be omitted to avoid duplicating the cassette."
                        ),
                    )
                )

        if requested_promoter_names and any(p in backbone_promoters for p in requested_promoter_names):
            include_promoter = False
            if insertion_mode == "full_cassette":
                insertion_mode = "cds_plus_tags"
            warnings.append(
                WarningMessage(
                    code="PROMOTER_ALREADY_IN_BACKBONE",
                    message="Requested promoter already appears to be present in the selected backbone; promoter will be omitted from the insert.",
                )
            )

        if inp.backbone.backbone_payload_markers:
            warnings.append(
                WarningMessage(
                    code="BACKBONE_PAYLOAD_MARKERS_DETECTED",
                    message=f"Backbone payload markers detected: {inp.backbone.backbone_payload_markers}",
                )
            )

        # Promoter
        promoter_included = False
        if include_promoter and promoters:
            p = promoters[0]
            cursor = self._append_feature(
                sequence_parts=sequence_parts,
                annotations=annotations,
                feature_order=feature_order,
                cursor=cursor,
                name=p.name,
                seq=p.sequence,
                feat_type="promoter",
                order_label=f"{p.name}_promoter",
            )
            promoter_included = True
        elif not promoters and insertion_mode == "full_cassette":
            warnings.append(
                WarningMessage(
                    code="NO_PROMOTER",
                    message="No promoter feature resolved; construct may not be expressed.",
                )
            )

        # Kozak
        if kozak_feats:
            if "bacterial" in host_system or "e. coli" in host_system:
                warnings.append(
                    WarningMessage(
                        code="KOZAK_IN_NON_MAMMALIAN_CONTEXT",
                        message="Kozak sequence present in a non-mammalian expression context.",
                    )
                )

            k = kozak_feats[0]
            cursor = self._append_feature(
                sequence_parts=sequence_parts,
                annotations=annotations,
                feature_order=feature_order,
                cursor=cursor,
                name="Kozak",
                seq=k.sequence,
                feat_type="misc_feature",
                order_label="Kozak",
            )

        # N-tags
        for tag in n_tags:
            cursor = self._append_feature(
                sequence_parts=sequence_parts,
                annotations=annotations,
                feature_order=feature_order,
                cursor=cursor,
                name=f"{tag.name}_tag",
                seq=tag.sequence,
                feat_type="misc_feature",
                order_label=f"{tag.name}_N_tag",
            )

        if unspecified_tags:
            warnings.append(
                WarningMessage(
                    code="TAG_POSITION_UNSPECIFIED",
                    message="One or more tag features have no position_hint. Unspecified tags will be placed after the CDS by default.",
                )
            )

        # CDS
        cds = self._clean_seq(inp.gene.cds_sequence or "")
        warnings.extend(self._cds_warnings(cds))

        if cds:
            gene_label = f"{inp.gene.gene_symbol}_CDS"
            cursor = self._append_feature(
                sequence_parts=sequence_parts,
                annotations=annotations,
                feature_order=feature_order,
                cursor=cursor,
                name=gene_label,
                seq=cds,
                feat_type="CDS",
                order_label=gene_label,
            )

        # Unspecified tags
        for tag in unspecified_tags:
            cursor = self._append_feature(
                sequence_parts=sequence_parts,
                annotations=annotations,
                feature_order=feature_order,
                cursor=cursor,
                name=f"{tag.name}_tag",
                seq=tag.sequence,
                feat_type="misc_feature",
                order_label=f"{tag.name}_tag",
            )

        # C-tags
        for tag in c_tags:
            cursor = self._append_feature(
                sequence_parts=sequence_parts,
                annotations=annotations,
                feature_order=feature_order,
                cursor=cursor,
                name=f"{tag.name}_tag",
                seq=tag.sequence,
                feat_type="misc_feature",
                order_label=f"{tag.name}_C_tag",
            )

        # Terminator / polyA
        terminator_included = False
        if include_terminator and terminators and insertion_mode == "full_cassette":
            t = terminators[0]
            cursor = self._append_feature(
                sequence_parts=sequence_parts,
                annotations=annotations,
                feature_order=feature_order,
                cursor=cursor,
                name=t.name,
                seq=t.sequence,
                feat_type="terminator",
                order_label=t.name,
            )
            terminator_included = True
        elif not terminators and "mammalian" in host_system and insertion_mode == "full_cassette":
            warnings.append(
                WarningMessage(
                    code="NO_TERMINATOR_OR_POLYA",
                    message="Mammalian expression construct has no terminator/polyA feature. The cassette may be incomplete.",
                )
            )

        construct_seq = "".join(sequence_parts)
        if not construct_seq:
            warnings.append(
                WarningMessage(
                    code="EMPTY_CONSTRUCT",
                    message="Construct sequence is empty.",
                )
            )

        out = ConstructOutput(
            feature_order=feature_order,
            construct_sequence=construct_seq,
            annotations=annotations,
            warnings=warnings,
            promoter_included=promoter_included,
            terminator_included=terminator_included,
            polyA_included=False,
            insertion_mode=insertion_mode,  # type: ignore[arg-type]
        )

        self._logger.info("[CONSTRUCT] OUTPUT %s", out.model_dump())
        return out