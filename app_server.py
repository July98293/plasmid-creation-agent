"""
Flask development server for the plasmid-creation frontend.

Serves:
  GET  /                          → app/index.html
  POST /api/intent                → IntentAgent (parse user request)
  POST /api/pipeline              → full pipeline, one blocking response
  POST /api/pipeline/start        → stores request, returns {job_id}
  GET  /api/pipeline/stream/<id>  → SSE stream with live element-by-element progress
"""
from __future__ import annotations

import json
import os
import queue as queue_module
import threading
import traceback
import uuid

from dotenv import load_dotenv
from flask import Flask, Response, jsonify, request, send_from_directory
from flask_cors import CORS

from pipeline_agent.intent_agent import IntentAgent
from pipeline_agent.models import IntentInput, OutputInput, PlasmidConstructionInput
from pipeline_agent.orchestrator import PipelineOrchestrator
from pipeline_agent.output_agent import OutputAgent
from pipeline_agent.plasmid_construction_agent import PlasmidConstructionAgent

load_dotenv()

app = Flask(__name__, static_folder="app", static_url_path="")
CORS(app)

# In-memory job store {job_id: user_request}
_jobs: dict[str, str] = {}


# ---------------------------------------------------------------------------
# Frontend
# ---------------------------------------------------------------------------

@app.route("/")
def index():
    return send_from_directory("app", "index.html")


# ---------------------------------------------------------------------------
# API: intent only (fast preview)
# ---------------------------------------------------------------------------

@app.route("/api/intent", methods=["POST"])
def api_intent():
    body = request.get_json(force=True) or {}
    user_request = (body.get("user_request") or "").strip()
    if not user_request:
        return jsonify({"error": "user_request is required"}), 400
    try:
        intent = IntentAgent().run(IntentInput(user_request=user_request))
        return jsonify(intent.model_dump())
    except Exception as exc:
        traceback.print_exc()
        return jsonify({"error": str(exc)}), 500


# ---------------------------------------------------------------------------
# API: full pipeline — single blocking response (kept for run_test.py)
# ---------------------------------------------------------------------------

@app.route("/api/pipeline", methods=["POST"])
def api_pipeline():
    body = request.get_json(force=True) or {}
    user_request = (body.get("user_request") or "").strip()
    if not user_request:
        return jsonify({"error": "user_request is required"}), 400
    try:
        result = PipelineOrchestrator().run(user_request)
        return jsonify({
            "intent": result.intent.model_dump(),
            "construction": result.construction.model_dump(),
            "export_output": {
                "gbk_text":  result.output.gbk_text,
                "file_name": result.output.file_name,
                "warnings":  result.output.frontend_warnings,
            },
        })
    except Exception as exc:
        traceback.print_exc()
        return jsonify({"error": str(exc)}), 500


# ---------------------------------------------------------------------------
# API: streaming pipeline — step 1: register the job
# ---------------------------------------------------------------------------

@app.route("/api/pipeline/start", methods=["POST"])
def api_pipeline_start():
    body = request.get_json(force=True) or {}
    user_request = (body.get("user_request") or "").strip()
    if not user_request:
        return jsonify({"error": "user_request is required"}), 400
    job_id = str(uuid.uuid4())
    _jobs[job_id] = user_request
    return jsonify({"job_id": job_id})


# ---------------------------------------------------------------------------
# API: streaming pipeline — step 2: SSE stream
# ---------------------------------------------------------------------------

def _sse(event: dict) -> str:
    return f"data: {json.dumps(event)}\n\n"


@app.route("/api/pipeline/stream/<job_id>")
def api_pipeline_stream(job_id: str):
    user_request = _jobs.pop(job_id, None)
    if not user_request:
        return jsonify({"error": "job not found"}), 404

    q: queue_module.Queue = queue_module.Queue()

    def worker():
        try:
            # ── 1. Intent ────────────────────────────────────────────────
            q.put({"type": "step", "label": "Interpreting request…"})
            intent = IntentAgent().run(IntentInput(user_request=user_request))
            q.put({"type": "intent", "data": intent.model_dump()})

            # ── 2. Construction (element-by-element via progress_cb) ─────
            q.put({"type": "step", "label": "Assembling plasmid elements…"})

            def progress_cb(event: dict):
                q.put({"type": "progress", **event})

            construction = PlasmidConstructionAgent().run(
                PlasmidConstructionInput(intent=intent),
                progress_cb=progress_cb,
            )
            q.put({"type": "construction", "data": construction.model_dump()})

            # ── 3. Export ────────────────────────────────────────────────
            q.put({"type": "step", "label": "Generating GenBank file…"})
            output = OutputAgent().run(
                OutputInput(construct_output=construction, name="plasmid")
            )
            q.put({
                "type": "complete",
                "data": {
                    "intent":       intent.model_dump(),
                    "construction": construction.model_dump(),
                    "export_output": {
                        "gbk_text":  output.gbk_text,
                        "file_name": output.file_name,
                        "warnings":  output.frontend_warnings,
                    },
                },
            })
        except Exception as exc:
            traceback.print_exc()
            q.put({"type": "error", "message": str(exc)})
        finally:
            q.put(None)  # sentinel

    threading.Thread(target=worker, daemon=True).start()

    def generate():
        while True:
            try:
                event = q.get(timeout=180)
            except queue_module.Empty:
                yield _sse({"type": "ping"})
                continue
            if event is None:
                break
            yield _sse(event)

    return Response(
        generate(),
        mimetype="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    port = int(os.getenv("PORT", 5000))
    print(f"\n  Plasmid UI →  http://localhost:{port}\n")
    app.run(host="0.0.0.0", port=port, debug=True, threaded=True)
