from __future__ import annotations

import re
from typing import Optional

from logging_utils import get_conversation_logger

from .models import IntentInput, IntentOutput


class IntentAgent:
    """
    Deterministic-ish intent parser.

    For now this is regex / rule based so the pipeline stays deterministic
    even without an external LLM, but the IO is Pydantic and stable.
    """

    def __init__(self) -> None:
        self._logger = get_conversation_logger()

    def _infer_species(self, text: str) -> str:
        t = text.lower()
        if any(k in t for k in ["hek", "293", "human", "mammal"]):
            return "homo_sapiens"
        if "mouse" in t or "murine" in t or "cho" in t:
            return "mus_musculus"
        if "yeast" in t or "saccharomyces" in t:
            return "saccharomyces_cerevisiae"
        if "ecoli" in t or "e. coli" in t or "e-coli" in t:
            return "escherichia_coli"
        return "homo_sapiens"

    def _infer_host(self, text: str) -> str:
        t = text.lower()
        if "hek" in t or "293" in t:
            return "HEK293"
        if "cho" in t:
            return "CHO"
        if "ecoli" in t or "e. coli" in t or "e-coli" in t:
            return "E. coli"
        if "yeast" in t or "saccharomyces" in t:
            return "yeast"
        return "mammalian (unspecified)"

    def _infer_promoter(self, text: str) -> Optional[str]:
        t = text.lower()
        if "cmv" in t:
            return "CMV"
        if "ef1a" in t or "ef1α" in t:
            return "EF1a"
        if "cag" in t:
            return "CAG"
        if "pgk" in t:
            return "PGK"
        if "t7" in t:
            return "T7"
        if "lac" in t:
            return "lac"
        return None

    def _infer_assembly(self, text: str) -> Optional[str]:
        t = text.lower()
        if "gibson" in t:
            return "Gibson"
        if "golden gate" in t or "goldengate" in t:
            return "GoldenGate"
        if "restriction" in t or "digest" in t or "ligat" in t:
            return "Restriction"
        if "synthesis" in t or "synthesise" in t:
            return "Synthesis"
        return None

    def _infer_tags(self, text: str) -> tuple[Optional[str], Optional[str]]:
        t = text.lower()
        n_tag: Optional[str] = None
        c_tag: Optional[str] = None
        if "flag tag" in t or "flag-tag" in t or "n-terminal flag" in t:
            if "c-terminal" in t:
                c_tag = "FLAG"
            else:
                n_tag = "FLAG"
        if "his tag" in t or "6xhis" in t or "his6" in t:
            c_tag = c_tag or "His6"
        return n_tag, c_tag

    def _infer_terminator(self, text: str) -> Optional[str]:
        t = text.lower()
        if "bgh" in t and "poly" in t:
            return "BGH_polyA"
        if "sv40" in t and "poly" in t:
            return "SV40_polyA"
        return None

    def _extract_gene_symbol(self, text: str) -> Optional[str]:
        # Very small heuristic: first all-caps token with 2–8 chars
        tokens = re.findall(r"\b[A-Z0-9]{2,8}\b", text)
        for tok in tokens:
            if tok in {"CMV", "T7"}:
                continue
            return tok
        return None

    def run(self, inp: IntentInput) -> IntentOutput:
        self._logger.info("[INTENT] INPUT %s", inp.user_request)

        text = inp.user_request

        gene_symbol = self._extract_gene_symbol(text) or "GENE"
        target_species = self._infer_species(text)
        expression_host = self._infer_host(text)
        promoter = self._infer_promoter(text)
        assembly_method = self._infer_assembly(text)
        n_tag, c_tag = self._infer_tags(text)
        terminator = self._infer_terminator(text)

        out = IntentOutput(
            gene_symbol=gene_symbol,
            target_species=target_species,
            expression_host=expression_host,
            promoter=promoter,
            assembly_method=assembly_method,
            n_terminal_tag=n_tag,
            c_terminal_tag=c_tag,
            terminator=terminator,
        )

        self._logger.info("[INTENT] OUTPUT %s", out.model_dump())
        return out

