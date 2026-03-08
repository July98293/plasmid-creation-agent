from __future__ import annotations

from logging_utils import get_conversation_logger

from .models import ExpressionInput, ExpressionOutput, WarningMessage


class ExpressionAgent:
    """
    Simple rule-based expression context checker.
    """

    def __init__(self) -> None:
        self._logger = get_conversation_logger()

    def _host_system(self, host: str) -> str:
        h = host.lower()
        if "hek" in h or "293" in h or "cho" in h or "mammal" in h:
            return "mammalian transient expression"
        if "e. coli" in h or "ecoli" in h or "e-coli" in h:
            return "bacterial expression (E. coli)"
        if "yeast" in h or "saccharomyces" in h:
            return "yeast expression"
        return "unspecified host system"

    def run(self, inp: ExpressionInput) -> ExpressionOutput:
        self._logger.info("[EXPRESSION] INPUT %s", inp.model_dump())

        host_system = self._host_system(inp.intent.expression_host)

        # Basic codon optimization heuristic: if species doesn't match host category.
        codon_needed = False
        species = inp.gene.organism
        if "homo_sapiens" in species and "bacterial" in host_system:
            codon_needed = True
        if "escherichia_coli" in species and "mammalian" in host_system:
            codon_needed = True

        warnings: list[WarningMessage] = []
        notes: list[str] = []

        promoter_names = [f.name for f in inp.features.features if f.type == "promoter"]
        if host_system.startswith("mammalian"):
            if any(p in {"CMV", "EF1a", "CAG", "PGK"} for p in promoter_names):
                notes.append("Promoter is suitable for mammalian expression.")
            else:
                warnings.append(
                    WarningMessage(
                        code="PROMOTER_HOST_MISMATCH",
                        message="Promoter may not be optimal for mammalian expression.",
                    )
                )

        if any(f.type == "kozak" for f in inp.features.features):
            notes.append("Kozak sequence present upstream of CDS.")
        else:
            warnings.append(
                WarningMessage(
                    code="KOZAK_MISSING",
                    message="No Kozak sequence detected; translation efficiency may be reduced.",
                )
            )

        if codon_needed:
            notes.append("Codon optimization recommended between gene species and host.")

        out = ExpressionOutput(
            host_system=host_system,
            codon_optimization_needed=codon_needed,
            compatibility_warnings=warnings,
            notes=notes,
        )

        self._logger.info("[EXPRESSION] OUTPUT %s", out.model_dump())
        return out

