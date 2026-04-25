from __future__ import annotations

import json
import sys
import traceback
from typing import Any

from dotenv import load_dotenv

from pipeline_agent.intent_agent import IntentAgent
from pipeline_agent.models import IntentInput, OutputInput, PlasmidConstructionInput
from pipeline_agent.orchestrator import PipelineOrchestrator
from pipeline_agent.output_agent import OutputAgent
from pipeline_agent.plasmid_construction_agent import PlasmidConstructionAgent

load_dotenv()

DEFAULT_REQUEST = "Create a fusion plasmid with Gene A TP53 and Gene B EGFR in HEK293 using CMV promoter"


def pretty(title: str, data: Any) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)
    if hasattr(data, "model_dump"):
        print(json.dumps(data.model_dump(), indent=2, ensure_ascii=False))
    else:
        print(json.dumps(data, indent=2, ensure_ascii=False, default=str))


def run_intent(request: str) -> None:
    intent = IntentAgent().run(IntentInput(user_request=request))
    pretty("IntentOutput", intent)


def run_construction(request: str) -> None:
    intent = IntentAgent().run(IntentInput(user_request=request))
    pretty("IntentOutput", intent)
    construction = PlasmidConstructionAgent().run(PlasmidConstructionInput(intent=intent))
    pretty("PlasmidConstructionOutput", construction)


def run_output(request: str) -> None:
    intent = IntentAgent().run(IntentInput(user_request=request))
    construction = PlasmidConstructionAgent().run(PlasmidConstructionInput(intent=intent))
    output = OutputAgent().run(OutputInput(construct_output=construction, name="test_construct"))
    pretty("Output metadata", {"file_name": output.file_name, "frontend_warnings": output.frontend_warnings})
    print("\n--- GBK preview (first 1000 chars) ---\n")
    print(output.gbk_text[:1000])


def run_pipeline(request: str) -> None:
    result = PipelineOrchestrator().run(request)
    pretty("IntentOutput", result.intent)
    pretty("PlasmidConstructionOutput", result.construction)
    pretty("Output metadata", {"file_name": result.output.file_name, "frontend_warnings": result.output.frontend_warnings})
    print("\n--- GBK preview (first 1000 chars) ---\n")
    print(result.output.gbk_text[:1000])


def main() -> None:
    args = sys.argv[1:]
    if not args:
        print("Usage: python run_test.py <agent> [request]")
        print("  intent        [request]")
        print("  construction  [request]")
        print("  output        [request]")
        print("  pipeline      [request]")
        sys.exit(1)

    agent_name = args[0].lower()
    request = args[1] if len(args) > 1 else DEFAULT_REQUEST

    print(f"\nTesting: {agent_name!r}")
    print(f"Request: {request!r}")

    try:
        if agent_name == "intent":
            run_intent(request)
        elif agent_name == "construction":
            run_construction(request)
        elif agent_name == "output":
            run_output(request)
        elif agent_name == "pipeline":
            run_pipeline(request)
        else:
            print(f"Unknown mode: {agent_name!r}")
            print("Available: intent, construction, output, pipeline")
            sys.exit(1)
    except Exception as exc:
        print(f"\nFAILED: {exc}")
        traceback.print_exc()


if __name__ == "__main__":
    main()
