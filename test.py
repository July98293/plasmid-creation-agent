from __future__ import annotations

import asyncio
import json
import traceback
from typing import Any

from plasmid_pipeline.intent_agent import IntentAgent
from plasmid_pipeline.gene_agent import GeneAgent
from plasmid_pipeline.orchestrator import PipelineOrchestrator, PipelineValidationError
from plasmid_pipeline.models import IntentInput, GeneInput


def pretty(title: str, data: Any) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)
    if hasattr(data, "model_dump"):
        print(json.dumps(data.model_dump(), indent=2, ensure_ascii=False))
    else:
        print(json.dumps(data, indent=2, ensure_ascii=False, default=str))


def short_seq(seq: str | None, max_len: int = 120) -> str | None:
    if not seq:
        return seq
    if len(seq) <= max_len:
        return seq
    return seq[:max_len] + f"... [len={len(seq)}]"


async def test_intent(user_request: str) -> None:
    print("\n\n##############################")
    print("TEST 1 — INTENT AGENT")
    print("##############################")

    agent = IntentAgent()

    try:
        result = agent.run(IntentInput(user_request=user_request))
        pretty("Intent output", result)
    except Exception as exc:
        print("IntentAgent FAILED")
        print(f"Error: {exc}")
        traceback.print_exc()


async def test_gene_from_request(user_request: str) -> None:
    print("\n\n##############################")
    print("TEST 2 — INTENT -> GENE")
    print("##############################")

    intent_agent = IntentAgent()
    gene_agent = GeneAgent()

    try:
        intent = intent_agent.run(IntentInput(user_request=user_request))
        pretty("Intent output", intent)

        gene_input = GeneInput(
            gene_symbol=intent.gene_symbol,
            target_species=intent.target_species,
        )
        pretty("Gene input", gene_input)

        gene = await gene_agent.run(gene_input)

        gene_dump = gene.model_dump()
        gene_dump["cds_sequence"] = short_seq(gene_dump.get("cds_sequence"))
        pretty("Gene output", gene_dump)

        if gene.warnings:
            print("\nWarnings:")
            for w in gene.warnings:
                print(f"- [{w.code}] {w.message}")

    except Exception as exc:
        print("GeneAgent FAILED")
        print(f"Error: {exc}")
        traceback.print_exc()


async def test_gene_direct(gene_symbol: str, target_species: str) -> None:
    print("\n\n##############################")
    print("TEST 3 — GENE DIRECT")
    print("##############################")

    gene_agent = GeneAgent()

    try:
        gene_input = GeneInput(
            gene_symbol=gene_symbol,
            target_species=target_species,
        )
        pretty("Gene input", gene_input)

        gene = await gene_agent.run(gene_input)

        gene_dump = gene.model_dump()
        gene_dump["cds_sequence"] = short_seq(gene_dump.get("cds_sequence"))
        pretty("Gene output", gene_dump)

        if gene.warnings:
            print("\nWarnings:")
            for w in gene.warnings:
                print(f"- [{w.code}] {w.message}")

    except Exception as exc:
        print("Gene direct FAILED")
        print(f"Error: {exc}")
        traceback.print_exc()


async def test_orchestrator(user_request: str) -> None:
    print("\n\n##############################")
    print("TEST 4 — FULL ORCHESTRATOR")
    print("##############################")

    orchestrator = PipelineOrchestrator()

    try:
        result = await orchestrator.run(user_request)

        pretty("Intent", result.intent)

        gene_dump = result.gene.model_dump()
        gene_dump["cds_sequence"] = short_seq(gene_dump.get("cds_sequence"))
        pretty("Gene", gene_dump)

        pretty("Features", result.features)
        pretty("Expression", result.expression)

        backbone_dump = result.backbone.model_dump()
        backbone_dump["backbone_sequence"] = short_seq(backbone_dump.get("backbone_sequence"))
        pretty("Backbone", backbone_dump)

        construct_dump = result.construct_output.model_dump()
        construct_dump["construct_sequence"] = short_seq(construct_dump.get("construct_sequence"))
        pretty("Construct", construct_dump)

        assembly_dump = result.assembly.model_dump()
        assembly_dump["assembled_sequence"] = short_seq(assembly_dump.get("assembled_sequence"))
        pretty("Assembly", assembly_dump)

        pretty("Export summary", result.export_output.summary_json)

    except PipelineValidationError as exc:
        print("PipelineValidationError")
        print(exc)
        traceback.print_exc()

    except Exception as exc:
        print("Orchestrator FAILED")
        print(f"Error: {exc}")
        traceback.print_exc()


async def main() -> None:
    # Cambia qui con le richieste che vuoi provare
    requests = [
        "Express TP53 in HEK293 cells using a CMV promoter assembled with Gibson",
        "Build a GFP plasmid for E. coli with a T7 promoter",
        "Express EGFR in human cells with EF1a promoter and C-terminal His tag",
    ]

    # 1) test solo intent
    await test_intent(requests[0])

    # 2) test intent -> gene
    await test_gene_from_request(requests[0])

    # 3) test gene diretto
    await test_gene_direct("TP53", "homo_sapiens")
    await test_gene_direct("EGFR", "homo_sapiens")
    await test_gene_direct("GFP", "homo_sapiens")  # utile per vedere casi che possono fallire

    # 4) test full pipeline
    for req in requests:
        await test_orchestrator(req)


if __name__ == "__main__":
    asyncio.run(main())