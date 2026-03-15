from __future__ import annotations

import asyncio
import json
import sys
import traceback
from typing import Any

from dotenv import load_dotenv
load_dotenv()

from plasmid_pipeline.intent_agent import IntentAgent
from plasmid_pipeline.gene_agent import GeneAgent
from plasmid_pipeline.backbone_agent import BackboneAgent
from plasmid_pipeline.feature_agent import FeatureAgent
from plasmid_pipeline.expression_agent import ExpressionAgent
from plasmid_pipeline.assembly_agent import AssemblyAgent
from plasmid_pipeline.models import (
    IntentInput, GeneInput, BackboneInput, FeatureInput,
    ExpressionInput, AssemblyInput,
)


DEFAULT_REQUEST = "Express TP53 in HEK293 cells using a CMV promoter assembled with Gibson"


def pretty(title: str, data: Any) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)
    if hasattr(data, "model_dump"):
        print(json.dumps(data.model_dump(), indent=2, ensure_ascii=False))
    else:
        print(json.dumps(data, indent=2, ensure_ascii=False, default=str))


async def run_intent(request: str) -> None:
    inp = IntentInput(user_request=request)
    pretty("INPUT — IntentInput", inp)

    agent = IntentAgent()
    result = agent.run(inp)
    pretty("OUTPUT — IntentOutput", result)


async def run_backbone(request: str) -> None:
    intent_inp = IntentInput(user_request=request)
    pretty("INPUT — IntentInput", intent_inp)

    intent_agent = IntentAgent()
    intent = intent_agent.run(intent_inp)
    pretty("OUTPUT — IntentOutput", intent)

    backbone_inp = BackboneInput(
        expression_host=intent.expression_host or "",
        backbone=intent.backbone,
        vector_type=None,
    )
    pretty("INPUT — BackboneInput", backbone_inp)

    backbone_agent = BackboneAgent()
    result = backbone_agent.run(backbone_inp)
    pretty("OUTPUT — BackboneOutput", result)


async def run_feature(request: str) -> None:
    intent_inp = IntentInput(user_request=request)
    pretty("INPUT — IntentInput", intent_inp)

    intent_agent = IntentAgent()
    intent = intent_agent.run(intent_inp)
    pretty("OUTPUT — IntentOutput", intent)

    feature_inp = FeatureInput(
        promoter=intent.promoter,
        n_terminal_tag=intent.n_terminal_tag,
        c_terminal_tag=intent.c_terminal_tag,
        extra_tags=intent.extra_tags,
        terminator=intent.terminator,
        polyA=intent.polyA,
        selection_marker=intent.selection_marker,
        origin_of_replication=intent.origin_of_replication,
        expression_host=intent.expression_host,
    )
    pretty("INPUT — FeatureInput", feature_inp)

    feature_agent = FeatureAgent()
    result = feature_agent.run(feature_inp)
    pretty("OUTPUT — FeatureOutput", result)


async def run_gene(request: str) -> None:
    intent_inp = IntentInput(user_request=request)
    pretty("INPUT — IntentInput", intent_inp)

    intent_agent = IntentAgent()
    intent = intent_agent.run(intent_inp)
    pretty("OUTPUT — IntentOutput", intent)

    gene_inp = GeneInput(
        gene_symbol=intent.gene_symbol or "",
        target_species=intent.target_species,
    )
    pretty("INPUT — GeneInput", gene_inp)

    gene_agent = GeneAgent()
    result = await gene_agent.run(gene_inp)
    pretty("OUTPUT — GeneOutput", result)


async def run_assembly(request: str) -> None:
    intent = IntentAgent().run(IntentInput(user_request=request))
    pretty("IntentOutput", intent)

    gene, backbone, feature = await asyncio.gather(
        GeneAgent().run(GeneInput(gene_symbol=intent.gene_symbol or "", target_species=intent.target_species)),
        asyncio.to_thread(BackboneAgent().run, BackboneInput(expression_host=intent.expression_host or "", backbone=intent.backbone)),
        asyncio.to_thread(FeatureAgent().run, FeatureInput(
            promoter=intent.promoter, n_terminal_tag=intent.n_terminal_tag,
            c_terminal_tag=intent.c_terminal_tag, extra_tags=intent.extra_tags,
            terminator=intent.terminator, polyA=intent.polyA,
            selection_marker=intent.selection_marker,
            origin_of_replication=intent.origin_of_replication,
            expression_host=intent.expression_host,
        )),
    )
    pretty("GeneOutput", gene)
    pretty("BackboneOutput", backbone)
    pretty("FeatureOutput", feature)

    expression = ExpressionAgent().run(ExpressionInput(intent=intent, gene=gene, features=feature))
    pretty("ExpressionOutput", expression)

    assembly = AssemblyAgent().run(AssemblyInput(
        backbone_sequence=backbone.backbone_sequence,
        backbone_name=backbone.backbone_name,
        backbone_genbank=backbone.backbone_genbank,
        cds_sequence=gene.cds_sequence or "",
        gene_symbol=gene.gene_symbol,
        features=feature.features,
        assembly_method=intent.assembly_method,
    ))
    pretty("AssemblyOutput", assembly)


async def main() -> None:
    args = sys.argv[1:]
    if not args:
        print("Usage: python run_test.py <agent> [request]")
        print("  intent   [request]")
        print("  backbone [request]")
        print("  gene     [request]")
        print("  feature  [request]")
        print("  assembly [request]")
        sys.exit(1)

    agent_name = args[0].lower()
    request = args[1] if len(args) > 1 else DEFAULT_REQUEST

    print(f"\nTesting agent: {agent_name!r}")
    print(f"Request: {request!r}")

    try:
        if agent_name == "intent":
            await run_intent(request)
        elif agent_name == "backbone":
            await run_backbone(request)
        elif agent_name == "gene":
            await run_gene(request)
        elif agent_name == "feature":
            await run_feature(request)
        elif agent_name == "assembly":
            await run_assembly(request)
        else:
            print(f"Unknown agent: {agent_name!r}")
            print("Available: intent, backbone, gene, feature, assembly")
            sys.exit(1)
    except Exception as exc:
        print(f"\nFAILED: {exc}")
        traceback.print_exc()


if __name__ == "__main__":
    asyncio.run(main())
