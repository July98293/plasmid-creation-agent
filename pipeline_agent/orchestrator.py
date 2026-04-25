from __future__ import annotations

from pydantic import BaseModel

from .intent_agent import IntentAgent
from .models import (
    IntentInput,
    IntentOutput,
    OutputInput,
    OutputOutput,
    PlasmidConstructionInput,
    PlasmidConstructionOutput,
)
from .output_agent import OutputAgent
from .plasmid_construction_agent import PlasmidConstructionAgent


class PipelineResult(BaseModel):
    intent: IntentOutput
    construction: PlasmidConstructionOutput
    output: OutputOutput


class PipelineOrchestrator:
    def __init__(self) -> None:
        self.intent_agent = IntentAgent()
        self.construction_agent = PlasmidConstructionAgent()
        self.output_agent = OutputAgent()

    def run(self, user_request: str) -> PipelineResult:
        intent = self.intent_agent.run(IntentInput(user_request=user_request))
        construction = self.construction_agent.run(PlasmidConstructionInput(intent=intent))
        output = self.output_agent.run(OutputInput(construct_output=construction, name="pipeline_plasmid"))
        return PipelineResult(intent=intent, construction=construction, output=output)
