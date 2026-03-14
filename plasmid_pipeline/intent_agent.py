from __future__ import annotations

import json
import os
from typing import List, Optional

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
    tags: List[TagSpec] = Field(default_factory=list)
    selection_marker: Optional[str] = None
    origin_of_replication: Optional[str] = None
    terminator: Optional[str] = None
    polyA: Optional[str] = None
    cloning_site: Optional[str] = None
    expression_host: Optional[str] = None
    best_expression_host: Optional[str] = None
    assembly_method: Optional[str] = None
    target_species: Optional[str] = None
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
            "If a field is not stated, return null.\n"
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
            "Output schema keys exactly:\n"
            "{\n"
            '  "gene": null,\n'
            '  "backbone": null,\n'
            '  "promoter": null,\n'
            '  "tags": [],\n'
            '  "selection_marker": null,\n'
            '  "origin_of_replication": null,\n'
            '  "terminator": null,\n'
            '  "polyA": null,\n'
            '  "cloning_site": null,\n'
            '  "expression_host": null,\n'
            '  "best_expression_host": null,\n'
            '  "assembly_method": null,\n'
            '  "target_species": null,\n'
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

    def run(self, inp: IntentInput) -> IntentOutput:
        self._logger.info("[INTENT] INPUT %s", inp.user_request)

        extracted = self._extract_with_llm(inp.user_request)

        n_tag, c_tag, tag_notes = self._pick_n_and_c_tags(extracted.tags)
        all_notes = list(extracted.notes) + tag_notes

        # Map extracted fields into your existing pipeline IntentOutput.
        # This assumes your current IntentOutput still expects these legacy keys.
        out = IntentOutput(
            gene_symbol=extracted.gene or "GENE",
            target_species=extracted.target_species or "unspecified",
            expression_host=extracted.expression_host or "unspecified",
            promoter=extracted.promoter,
            assembly_method=extracted.assembly_method,
            n_terminal_tag=n_tag,
            c_terminal_tag=c_tag,
            terminator=extracted.terminator or extracted.polyA,
        )

        self._logger.info(
            "[INTENT] OUTPUT %s | extra_fields=%s",
            out.model_dump(),
            {
                "backbone": extracted.backbone,
                "selection_marker": extracted.selection_marker,
                "origin_of_replication": extracted.origin_of_replication,
                "polyA": extracted.polyA,
                "cloning_site": extracted.cloning_site,
                "best_expression_host": extracted.best_expression_host,
                "notes": all_notes,
            },
        )
        return out