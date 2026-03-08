from __future__ import annotations

from typing import List, Tuple

from logging_utils import get_conversation_logger

from .models import AssemblyInput, AssemblyOutput, AssemblyFragment, WarningMessage


class AssemblyAgent:
    """
    Assembly planner with a simple replacement model.

    Current behavior:
    - Gibson: replace a short region around a chosen insertion site and add
      homology arms to the insert fragment.
    - Golden Gate / Restriction / Synthesis: still simplified, but reported
      explicitly as such.

    This is still not a full plasmid CAD engine, but it is more honest than a
    pure midpoint insertion.
    """

    def __init__(self) -> None:
        self._logger = get_conversation_logger()

    def _clean_seq(self, seq: str) -> str:
        return "".join(ch for ch in (seq or "").upper() if ch in {"A", "T", "G", "C", "N"})

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

    def _choose_replacement_window(self, backbone: str, insert: str) -> Tuple[int, int, int]:
        """
        Returns:
        - cut_start (0-based inclusive)
        - cut_end   (0-based exclusive)
        - homology_len

        Strategy:
        - choose a central region
        - for larger backbones, replace a modest window rather than inserting blindly
        """
        blen = len(backbone)

        if blen == 0:
            return 0, 0, 0

        # Homology arm length
        if blen >= 4000:
            homology_len = 30
        elif blen >= 1500:
            homology_len = 25
        else:
            homology_len = 20

        center = blen // 2

        # Replace a modest chunk so final construct is "region replacement"
        # instead of "just shove insert into middle".
        replace_len = min(max(40, len(insert) // 20), 250)

        cut_start = max(homology_len, center - (replace_len // 2))
        cut_end = min(blen - homology_len, cut_start + replace_len)

        if cut_end <= cut_start:
            cut_start = center
            cut_end = center

        return cut_start, cut_end, homology_len

    def _local_duplication_warning(
        self,
        backbone: str,
        insert: str,
        cut_start: int,
        cut_end: int,
    ) -> List[WarningMessage]:
        warnings: List[WarningMessage] = []

        if not backbone or not insert:
            return warnings

        if insert in backbone:
            warnings.append(
                WarningMessage(
                    code="INSERT_ALREADY_PRESENT_IN_BACKBONE",
                    message=(
                        "The full insert sequence already appears inside the selected backbone. "
                        "This suggests the backbone may already contain related cargo."
                    ),
                )
            )

        # Check whether the first ~80 bp of insert already appear near the chosen region
        probe = insert[:80]
        local_start = max(0, cut_start - 300)
        local_end = min(len(backbone), cut_end + 300)
        local_context = backbone[local_start:local_end]

        if probe and len(probe) >= 40 and probe in local_context:
            warnings.append(
                WarningMessage(
                    code="LOCAL_JUNCTION_DUPLICATION",
                    message=(
                        "The 5' region of the insert already appears near the chosen assembly locus. "
                        "This may indicate promoter/cassette duplication."
                    ),
                )
            )

        return warnings

    def run(self, inp: AssemblyInput) -> AssemblyOutput:
        self._logger.info("[ASSEMBLY] INPUT %s", inp.model_dump())

        method = self._pick_method(inp.assembly_preference or "")

        insert = self._clean_seq(inp.construct_sequence)
        backbone = self._clean_seq(inp.backbone_sequence)

        fragments: List[AssemblyFragment] = []
        warnings: List[WarningMessage] = []

        if not backbone:
            warnings.append(
                WarningMessage(
                    code="EMPTY_BACKBONE",
                    message="Backbone sequence is empty; assembly output will be insert-only or empty.",
                )
            )

        if not insert:
            warnings.append(
                WarningMessage(
                    code="EMPTY_INSERT",
                    message="Construct sequence is empty; assembled plasmid equals backbone only.",
                )
            )

        if backbone and len(backbone) < 1000:
            warnings.append(
                WarningMessage(
                    code="SHORT_BACKBONE",
                    message="Backbone is shorter than 1000 bp; assembly context may be unrealistic.",
                )
            )

        if method == "Gibson" and backbone:
            cut_start, cut_end, hom_len = self._choose_replacement_window(backbone, insert)

            left_homology = backbone[max(0, cut_start - hom_len):cut_start]
            right_homology = backbone[cut_end:min(len(backbone), cut_end + hom_len)]

            if len(left_homology) < 15 or len(right_homology) < 15:
                warnings.append(
                    WarningMessage(
                        code="SHORT_HOMOLOGY",
                        message="Gibson homology arms shorter than 15 bp; assembly efficiency may be low.",
                    )
                )

            warnings.extend(
                self._local_duplication_warning(
                    backbone=backbone,
                    insert=insert,
                    cut_start=cut_start,
                    cut_end=cut_end,
                )
            )

            # Region replacement model
            backbone_left = backbone[:cut_start]
            backbone_right = backbone[cut_end:]
            assembled_sequence = backbone_left + insert + backbone_right

            fragments.append(
                AssemblyFragment(
                    name="backbone_left_arm",
                    sequence=backbone_left,
                    left_homology=None,
                    right_homology=left_homology if left_homology else None,
                )
            )

            fragments.append(
                AssemblyFragment(
                    name="insert_fragment",
                    sequence=insert,
                    left_homology=left_homology if left_homology else None,
                    right_homology=right_homology if right_homology else None,
                )
            )

            fragments.append(
                AssemblyFragment(
                    name="backbone_right_arm",
                    sequence=backbone_right,
                    left_homology=right_homology if right_homology else None,
                    right_homology=None,
                )
            )

            primer_requirements: List[str] = [
                f"Forward primer with ~{len(left_homology)} bp 5' homology to the left backbone arm.",
                f"Reverse primer with ~{len(right_homology)} bp 3' homology to the right backbone arm.",
                "Linearize or PCR-amplify the backbone to remove the replaced region before Gibson assembly.",
            ]

            junction_offset = len(backbone_left) + 1

        elif method == "Golden Gate":
            assembled_sequence = backbone + insert
            fragments.append(
                AssemblyFragment(
                    name="backbone_fragment",
                    sequence=backbone,
                    left_homology=None,
                    right_homology=None,
                )
            )
            fragments.append(
                AssemblyFragment(
                    name="insert_fragment",
                    sequence=insert,
                    left_homology=None,
                    right_homology=None,
                )
            )
            primer_requirements = [
                "Primers or fragments with flanking type IIS restriction sites.",
                "Ensure correct overhang design and absence of forbidden internal type IIS sites.",
            ]
            warnings.append(
                WarningMessage(
                    code="SIMPLIFIED_GOLDEN_GATE_MODEL",
                    message="Golden Gate assembly is currently modeled as simplified concatenation, not full overhang-aware ligation.",
                )
            )
            junction_offset = (len(backbone) + 1) if backbone else 1

        elif method == "Restriction":
            assembled_sequence = backbone + insert
            fragments.append(
                AssemblyFragment(
                    name="backbone_fragment",
                    sequence=backbone,
                    left_homology=None,
                    right_homology=None,
                )
            )
            fragments.append(
                AssemblyFragment(
                    name="insert_fragment",
                    sequence=insert,
                    left_homology=None,
                    right_homology=None,
                )
            )
            primer_requirements = [
                "Insert and backbone must have compatible restriction sites.",
                "Verify that no forbidden internal restriction sites disrupt cloning.",
            ]
            warnings.append(
                WarningMessage(
                    code="SIMPLIFIED_RESTRICTION_MODEL",
                    message="Restriction cloning is currently modeled as simplified concatenation, not enzyme-site-aware ligation.",
                )
            )
            junction_offset = (len(backbone) + 1) if backbone else 1

        else:  # Synthesis or fallback
            assembled_sequence = backbone + insert
            fragments.append(
                AssemblyFragment(
                    name="full_backbone_fragment",
                    sequence=backbone,
                    left_homology=None,
                    right_homology=None,
                )
            )
            fragments.append(
                AssemblyFragment(
                    name="synthetic_insert_fragment",
                    sequence=insert,
                    left_homology=None,
                    right_homology=None,
                )
            )
            primer_requirements = [
                "No PCR primers strictly required if the full construct is synthesized by a provider.",
            ]
            warnings.append(
                WarningMessage(
                    code="SIMPLIFIED_SYNTHESIS_MODEL",
                    message="Synthesis mode currently returns simple backbone+insert composition for planning purposes.",
                )
            )
            junction_offset = (len(backbone) + 1) if backbone else 1

        out = AssemblyOutput(
            assembly_method=method,
            fragments=fragments,
            primer_requirements=primer_requirements,
            assembled_sequence=assembled_sequence,
            junction_offset=junction_offset,
            warnings=warnings,
        )

        self._logger.info("[ASSEMBLY] OUTPUT %s", out.model_dump())
        return out