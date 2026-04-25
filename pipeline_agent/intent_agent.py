from __future__ import annotations

import json
import os
import re
from typing import List, Optional, cast

from openai import OpenAI
from pydantic import BaseModel, Field, ValidationError

from logging_utils import get_conversation_logger
from .models import IntentInput, IntentOutput


class TagSpec(BaseModel):
    name: str
    position: Optional[str] = None  # "N-terminal" | "C-terminal" | None


class IntentExtraction(BaseModel):
    gene: Optional[str] = None
    backbone: Optional[str] = None
    promoter: Optional[str] = None
    promoter2: Optional[str] = None
    rbs: Optional[str] = None
    rbs2: Optional[str] = None
    tags: List[TagSpec] = Field(default_factory=list)
    selection_marker: Optional[str] = None
    origin_of_replication: Optional[str] = None
    terminator: Optional[str] = None
    polyA: Optional[str] = None
    cloning_site: Optional[str] = None
    regulatory_element: Optional[str] = None
    construct_type: Optional[str] = None  # standard | fusion | multi_cassette_fusion | synthetic_circuit
    gene_a: Optional[str] = None
    gene_b: Optional[str] = None
    linker_sequence: Optional[str] = None
    expression_host: Optional[str] = None
    best_expression_host: Optional[str] = None
    assembly_method: Optional[str] = None
    target_species: Optional[str] = None
    # Fusion chain (for multi-domain fusions beyond gene_a/gene_b)
    cleavage_sites: List[str] = Field(default_factory=list)
    # Synthetic circuit fields
    circuit_topology: Optional[str] = None        # e.g. "toggle switch", "AND gate"
    circuit_logic: Optional[str] = None           # e.g. "NOT", "AND"
    input_signals: List[str] = Field(default_factory=list)
    output_reporters: List[str] = Field(default_factory=list)
    circuit_regulators: List[str] = Field(default_factory=list)
    notes: List[str] = Field(default_factory=list)


class IntentAgent:
    """
    LLM-based intent parser.

    Reads the user's plasmid-design request and extracts structured design fields.
    Uses OpenAI Responses API, then validates with Pydantic before converting to
    the project's IntentOutput model.
    """

    def __init__(
        self,
        *,
        openai_api_key: Optional[str] = None,
        model: str = "gpt-5",
    ) -> None:
        self._logger = get_conversation_logger()
        api_key = openai_api_key or os.getenv("OPENAI_API_KEY")
        if not api_key:
            raise ValueError("OPENAI_API_KEY is not set.")
        self._client = OpenAI(api_key=api_key)
        self._model = model

    def _build_system_prompt(self) -> str:
        return (
            "You extract plasmid design intent from a user's free-text request.\n"
            "Return JSON only.\n"
            "Do not include markdown.\n"
            "Do not invent highly specific biological parts unless strongly implied.\n"
            "Extract from the conversation all the features that are not mentioned but required for the payload.\n"
            "If a field is not specified and you have to fill it in, add '(suggested)' at the end of the field\n"
            "tags must be a list of objects with keys: name, position.\n"
            "Valid position values are 'N-terminal', 'C-terminal', or null.\n"
            "expression_host = what the user explicitly requested.\n"
            "best_expression_host = your biological recommendation based on the request.\n"
            "origin_of_replication means plasmid origin / ori.\n"
            "polyA means polyadenylation signal.\n"
            "cloning_site means MCS / cloning site / insertion site if mentioned.\n"
            "target_species means the species of the gene/cargo if stated or strongly implied.\n"
            "assembly_method should capture things like Gibson, GoldenGate, Restriction, Synthesis.\n"
            "Keep standard gene symbols exactly as written when possible.\n"
            "\n"
            "construct_type rules:\n"
            "  - 'standard'             : single gene of interest\n"
            "  - 'fusion'               : two genes joined by a linker (gene_a + linker + gene_b)\n"
            "  - 'multi_cassette_fusion': two separate expression cassettes on the same backbone\n"
            "  - 'synthetic_circuit'    : genetic circuit with logic (toggle switches, oscillators,\n"
            "                             AND/NOT gates, etc.). Populate circuit_* fields.\n"
            "\n"
            "For synthetic_circuit, fill:\n"
            "  circuit_topology  → e.g. 'toggle switch', 'bistable', 'AND gate', 'oscillator'\n"
            "  circuit_logic     → e.g. 'AND', 'NOT', 'OR', 'NAND'\n"
            "  input_signals     → list of inducers/signals, e.g. ['aTc', 'IPTG']\n"
            "  output_reporters  → list of reporter genes, e.g. ['GFP', 'mCherry']\n"
            "  circuit_regulators→ list of regulatory proteins the circuit encodes, e.g. ['TetR', 'LacI']\n"
            "\n"
            "For fusion constructs with cleavage sites, list them in cleavage_sites.\n"
            "\n"
            "Output schema keys exactly:\n"
            "{\n"
            '  "gene": null,\n'
            '  "backbone": null,\n'
            '  "promoter": null,\n'
            '  "promoter2": null,\n'
            '  "rbs": null,\n'
            '  "rbs2": null,\n'
            '  "tags": [],\n'
            '  "selection_marker": null,\n'
            '  "origin_of_replication": null,\n'
            '  "terminator": null,\n'
            '  "polyA": null,\n'
            '  "cloning_site": null,\n'
            '  "regulatory_element": null,\n'
            '  "construct_type": "standard",\n'
            '  "gene_a": null,\n'
            '  "gene_b": null,\n'
            '  "linker_sequence": null,\n'
            '  "expression_host": null,\n'
            '  "best_expression_host": null,\n'
            '  "assembly_method": null,\n'
            '  "target_species": null,\n'
            '  "cleavage_sites": [],\n'
            '  "circuit_topology": null,\n'
            '  "circuit_logic": null,\n'
            '  "input_signals": [],\n'
            '  "output_reporters": [],\n'
            '  "circuit_regulators": [],\n'
            '  "notes": []\n'
            "}"
        )

    def _extract_with_llm(self, user_text: str) -> IntentExtraction:
        response = self._client.responses.create(
            model=self._model,
            input=[
                {
                    "role": "system",
                    "content": self._build_system_prompt(),
                },
                {
                    "role": "user",
                    "content": user_text,
                },
            ],
        )

        text = getattr(response, "output_text", None)
        if not text:
            raise ValueError("OpenAI returned no output_text for intent extraction.")

        try:
            data = json.loads(text)
        except json.JSONDecodeError as e:
            raise ValueError(f"Intent extraction did not return valid JSON. Raw output: {text}") from e

        try:
            return IntentExtraction.model_validate(data)
        except ValidationError as e:
            raise ValueError(f"Intent extraction JSON failed schema validation: {e}") from e

    def _pick_n_and_c_tags(
        self,
        tags: List[TagSpec],
    ) -> tuple[Optional[str], Optional[str], List[str]]:
        n_tag: Optional[str] = None
        c_tag: Optional[str] = None
        notes: List[str] = []

        for tag in tags:
            if tag.position == "N-terminal" and n_tag is None:
                n_tag = tag.name
            elif tag.position == "C-terminal" and c_tag is None:
                c_tag = tag.name
            else:
                notes.append(
                    f"Unplaced or additional tag extracted: {tag.name} (position={tag.position})."
                )

        return n_tag, c_tag, notes

    @staticmethod
    def _clean(value: Optional[str]) -> tuple[Optional[str], bool]:
        """Strip (suggested) suffix. Returns (cleaned_value, was_suggested)."""
        if value is None:
            return None, False
        was_suggested = bool(re.search(r"\(suggested\)", value, re.IGNORECASE))
        cleaned = re.sub(r"\s*\(suggested\)\s*$", "", value, flags=re.IGNORECASE).strip()
        return (cleaned or None), was_suggested

    def run(self, inp: IntentInput) -> IntentOutput:
        self._logger.info("[INTENT] INPUT %s", inp.user_request)

        extracted = self._extract_with_llm(inp.user_request)

        n_tag, c_tag, tag_notes = self._pick_n_and_c_tags(extracted.tags)
        all_notes = list(extracted.notes) + tag_notes

        suggested_fields: List[str] = []

        def cs(field_name: str, value: Optional[str]) -> Optional[str]:
            cleaned, was_suggested = self._clean(value)
            if was_suggested:
                suggested_fields.append(field_name)
            return cleaned

        # Also clean tag names (from _pick_n_and_c_tags, which returned raw values)
        n_tag_clean, n_suggested = self._clean(n_tag)
        c_tag_clean, c_suggested = self._clean(c_tag)
        if n_suggested:
            suggested_fields.append("n_terminal_tag")
        if c_suggested:
            suggested_fields.append("c_terminal_tag")

        _VALID_CONSTRUCT_TYPES = {"standard", "fusion", "multi_cassette_fusion", "synthetic_circuit"}
        construct_type_raw = (cs("construct_type", extracted.construct_type) or "standard").lower()
        if construct_type_raw not in _VALID_CONSTRUCT_TYPES:
            all_notes.append(
                f"Unrecognized construct_type '{construct_type_raw}', defaulted to 'standard'."
            )
            construct_type_raw = "standard"

        out = IntentOutput(
            gene_symbol=cs("gene_symbol", extracted.gene) or "GENE",
            target_species=cs("target_species", extracted.target_species) or "unspecified",
            expression_host=cs("expression_host", extracted.expression_host) or "unspecified",
            best_expression_host=cs("best_expression_host", extracted.best_expression_host),
            backbone=cs("backbone", extracted.backbone),
            promoter=cs("promoter", extracted.promoter),
            promoter2=cs("promoter2", extracted.promoter2),
            rbs=cs("rbs", extracted.rbs),
            rbs2=cs("rbs2", extracted.rbs2),
            assembly_method=cs("assembly_method", extracted.assembly_method),
            n_terminal_tag=n_tag_clean,
            c_terminal_tag=c_tag_clean,
            selection_marker=cs("selection_marker", extracted.selection_marker),
            origin_of_replication=cs("origin_of_replication", extracted.origin_of_replication),
            terminator=cs("terminator", extracted.terminator),
            polyA=cs("polyA", extracted.polyA),
            cloning_site=cs("cloning_site", extracted.cloning_site),
            regulatory_element=cs("regulatory_element", extracted.regulatory_element),
            construct_type=cast(str, construct_type_raw),  # validated by IntentOutput
            gene_a=cs("gene_a", extracted.gene_a),
            gene_b=cs("gene_b", extracted.gene_b),
            linker_sequence=cs("linker_sequence", extracted.linker_sequence),
            # Fusion chain cleavage sites
            cleavage_sites=extracted.cleavage_sites,
            # Synthetic circuit fields
            circuit_topology=cs("circuit_topology", extracted.circuit_topology),
            circuit_logic=cs("circuit_logic", extracted.circuit_logic),
            input_signals=extracted.input_signals,
            output_reporters=extracted.output_reporters,
            circuit_regulators=extracted.circuit_regulators,
            notes=all_notes,
            suggested_fields=suggested_fields,
        )

        self._logger.info("[INTENT] OUTPUT %s", out.model_dump())
        return out
