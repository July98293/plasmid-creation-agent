from __future__ import annotations

import asyncio

from logging_utils import get_conversation_logger

from .assembly_agent import AssemblyAgent
from .backbone_agent import BackboneAgent
from .construct_agent import ConstructAgent
from .expression_agent import ExpressionAgent
from .export_agent import ExportAgent
from .gene_agent import GeneAgent
from .intent_agent import IntentAgent
from .feature_agent import FeatureAgent
from .models import (
    IntentInput,
    GeneInput,
    FeatureInput,
    ExpressionInput,
    ConstructInput,
    AssemblyInput,
    ExportInput,
    BackboneInput,
    PipelineResult,
)


class PipelineValidationError(Exception):
    pass


class PipelineOrchestrator:
    """
    Deterministic pipeline:

    Intent -> Gene -> Feature -> Expression -> Backbone -> Construct -> Assembly -> Export
    """

    def __init__(self) -> None:
        self._logger = get_conversation_logger()

        self.intent_agent = IntentAgent()
        self.gene_agent = GeneAgent()
        self.feature_agent = FeatureAgent()
        self.expression_agent = ExpressionAgent()
        self.backbone_agent = BackboneAgent()
        self.construct_agent = ConstructAgent()
        self.assembly_agent = AssemblyAgent()
        self.export_agent = ExportAgent()

    def _require(self, condition: bool, message: str) -> None:
        if not condition:
            raise PipelineValidationError(message)

    async def run(self, user_request: str) -> PipelineResult:
        self._logger.info("[PIPELINE] START user_request=%s", user_request)

        intent = self.intent_agent.run(IntentInput(user_request=user_request))

        gene = await self.gene_agent.run(
            GeneInput(
                gene_symbol=intent.gene_symbol,
                target_species=intent.target_species,
            )
        )

        if not gene.cds_sequence:
            raise PipelineValidationError(
                f"No CDS sequence found for {intent.gene_symbol}; "
                f"gene_id={gene.gene_id}, transcript_id={gene.transcript_id}, "
                f"warnings={[w.model_dump() for w in gene.warnings]}"
            )

        feature = self.feature_agent.run(
            FeatureInput(
                promoter=intent.promoter,
                n_terminal_tag=intent.n_terminal_tag,
                c_terminal_tag=intent.c_terminal_tag,
                terminator=intent.terminator,
            )
        )

        expression = self.expression_agent.run(
            ExpressionInput(intent=intent, gene=gene, features=feature)
        )

        backbone = self.backbone_agent.run(
            BackboneInput(
                expression_host=intent.expression_host,
                promoter=intent.promoter,
                vector_type=None,
            )
        )

        construct = self.construct_agent.run(
            ConstructInput(
                gene=gene,
                features=feature,
                expression=expression,
                backbone=backbone,
            )
        )

        self._require(bool(construct.construct_sequence), "Construct sequence is empty.")

        if backbone.is_loaded_vector and construct.promoter_included:
            raise PipelineValidationError(
                "Selected backbone looks like a loaded expression vector, but the construct still includes a promoter. "
                "This would likely create a duplicate expression cassette."
            )

        assembly = self.assembly_agent.run(
            AssemblyInput(
                construct_sequence=construct.construct_sequence,
                backbone_sequence=backbone.backbone_sequence,
                assembly_preference=intent.assembly_method,
            )
        )

        export = self.export_agent.run(
            ExportInput(
                intent=intent,
                gene=gene,
                features=feature,
                expression=expression,
                construct_output=construct,
                backbone=backbone,
                assembly=assembly,
            )
        )

        result = PipelineResult(
            intent=intent,
            gene=gene,
            features=feature,
            expression=expression,
            backbone=backbone,
            construct_output=construct,
            assembly=assembly,
            export_output=export,
        )

        self._logger.info("[PIPELINE] END")
        return result


async def example_run() -> None:
    orchestrator = PipelineOrchestrator()
    req = "Express TP53 in HEK293 cells using a CMV promoter assembled with Gibson"
    result = await orchestrator.run(req)
    print("=== Summary ===")
    print(result.export_output.summary_json)
    print("\n=== Plasmid map data ===")
    print(result.export_output.plasmid_map_data)




if __name__ == "__main__":
    asyncio.run(example_run())