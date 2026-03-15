from __future__ import annotations

import asyncio
from pathlib import Path
from typing import Any, Dict

from flask import Flask, jsonify, request, send_from_directory
from flask_cors import CORS

from .intent_agent import IntentAgent
from .models import IntentInput
from .orchestrator import PipelineOrchestrator, PipelineValidationError


app = Flask(__name__, static_folder=None)
CORS(app)

_orchestrator = PipelineOrchestrator()
_intent_agent = IntentAgent()


def _run_pipeline(user_request: str) -> Dict[str, Any]:
    """
    Synchronous helper that runs the asyncio-based pipeline and returns a
    JSON-serialisable dict focused on export-layer fields.
    """
    result = asyncio.run(_orchestrator.run(user_request))

    export = result.export_output

    return {
        "intent": result.intent.model_dump(),
        "gene": result.gene.model_dump(),
        "features": result.features.model_dump(),
        "expression": result.expression.model_dump(),
        "backbone": result.backbone.model_dump(),
        "assembly": result.assembly.model_dump(),
        "export_output": {
            "summary_json": export.summary_json,
            "genbank_record": export.genbank_record,
            "plasmid_map_data": export.plasmid_map_data,
            "primer_list": export.primer_list,
            "bill_of_materials": export.bill_of_materials,
        },
    }


@app.route("/", methods=["GET"])
def index() -> Any:
    """
    Serve the single-page UI from the existing app/index.html file.
    """
    app_dir = Path(__file__).resolve().parent.parent / "app"
    return send_from_directory(app_dir, "index.html")


@app.route("/api/intent", methods=["POST"])
def run_intent() -> Any:
    """
    Run only the IntentAgent and return the structured intent with suggested_fields.

    Request JSON: { "user_request": "..." }
    """
    payload = request.get_json(silent=True) or {}
    user_request = str(payload.get("user_request") or "").strip()
    if not user_request:
        return jsonify({"error": "Missing 'user_request' in JSON body."}), 400

    try:
        intent = _intent_agent.run(IntentInput(user_request=user_request))
    except Exception as exc:
        return jsonify({"error": f"Intent extraction failed: {exc}"}), 500

    return jsonify(intent.model_dump())


@app.route("/api/pipeline", methods=["POST"])
def run_pipeline() -> Any:
    """
    Run the plasmid design pipeline for a free-text user_request.

    Request JSON:
    {
      "user_request": "Express TP53 in HEK293 cells using a CMV promoter assembled with Gibson"
    }
    """
    payload = request.get_json(silent=True) or {}
    user_request = str(payload.get("user_request") or "").strip()
    if not user_request:
        return jsonify({"error": "Missing 'user_request' in JSON body."}), 400

    try:
        data = _run_pipeline(user_request)
    except PipelineValidationError as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:  # noqa: BLE001
        return jsonify({"error": f"Unexpected server error: {exc}"}), 500

    return jsonify(data)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5167, debug=True)