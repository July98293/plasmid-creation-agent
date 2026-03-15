# Plasmid Creation Agent

A multi-agent pipeline that turns a free-text plasmid design request into a fully sequenced, annotated construct — complete with primer designs, a GenBank record, and a plasmid map.

---

## How it works

```
User request (natural language)
        │
        ▼
  ┌─────────────┐
  │ IntentAgent │  LLM (GPT) — extracts gene, host, promoter, backbone,
  └──────┬──────┘  assembly method, tags, selection marker, ori, polyA
         │
         ▼  (concurrent)
  ┌──────┴──────────────────────────────────────┐
  │  GeneAgent   BackboneAgent   FeatureAgent   │
  │  (async)     (thread)        (thread)       │
  │  Ensembl/    NCBI nuccore    NCBI + Addgene │
  │  NCBI        GenBank fetch   GenBank fetch  │
  └──────┬──────────────────────────────────────┘
         │
         ▼
  ┌──────────────────┐
  │ ExpressionAgent  │  rule-based host/codon compatibility check
  └──────┬───────────┘
         │
         ▼
  ┌───────────────┐
  │ ConstructAgent│  assembles cassette: promoter→kozak→tags→CDS→terminator
  └──────┬────────┘
         │
         ▼
  ┌───────────────┐
  │ AssemblyAgent │  designs Gibson / Golden Gate / restriction fragments + primers
  └──────┬────────┘
         │
         ▼
  ┌─────────────┐
  │ ExportAgent │  outputs GenBank record, plasmid map, primer list, BOM
  └─────────────┘
```

Gene, Backbone, and Feature agents run **in parallel** via `asyncio.gather` + `asyncio.to_thread` after Intent finishes.

---

## Agents

| Agent | Input | What it does | External APIs |
|---|---|---|---|
| **IntentAgent** | free text | LLM extraction of all design fields | OpenAI |
| **GeneAgent** | gene symbol, species | fetches CDS from Ensembl REST → NCBI fallback | Ensembl, NCBI |
| **BackboneAgent** | expression host, backbone name | searches NCBI nuccore for an empty scaffold, scores by size and name match | NCBI |
| **FeatureAgent** | promoter, polyA, terminator, selection marker, ori, tags | fetches sequences from NCBI (primary) then Addgene (fallback); Kozak is hardcoded | NCBI, Addgene |
| **ExpressionAgent** | intent + gene + features | checks host/codon compatibility, validates promoter suitability | — |
| **ConstructAgent** | gene + features + expression + backbone | builds linear cassette in biological order | — |
| **AssemblyAgent** | construct sequence + backbone sequence | designs assembly fragments, homology arms, primer requirements | — |
| **ExportAgent** | all outputs | produces JSON summary, GenBank record, plasmid map data, primer list, bill of materials | — |

---

## Setup

### Requirements

- Python 3.10+
- A `.env` file in the project root

### Install

```bash
python3 -m venv venv
source venv/bin/activate      # Windows: venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Environment variables

Create a `.env` file:

```env
OPENAI_API_KEY=sk-...           # required — used by IntentAgent
ADDGENE_API_KEY=...             # required — used by FeatureAgent fallback
NCBI_API_KEY=...                # optional — raises NCBI rate limit from 3 to 10 req/s
NCBI_EMAIL=you@example.com      # optional — NCBI best practice
AGENT_LOG_PATH=agent_conversation.log  # optional — defaults to cwd
```

---

## Testing individual agents

Use `run_test.py` to test any single agent in isolation. Each command runs **IntentAgent first** to extract structured inputs from your request, then pipes them into the target agent.

```bash
python run_test.py <agent> ["your request"]
```

The second argument is optional — it defaults to:
> `"Express TP53 in HEK293 cells using a CMV promoter assembled with Gibson"`

### Intent agent

Parses a free-text request and shows all extracted fields.

```bash
python run_test.py intent
python run_test.py intent "Clone EGFR into CHO cells with EF1a promoter"
```

### Gene agent

Runs Intent then fetches the CDS sequence for the extracted gene.

```bash
python run_test.py gene
python run_test.py gene "Express BRCA1 in HEK293 using pcDNA3.1"
```

### Backbone agent

Runs Intent then searches NCBI for an appropriate empty vector scaffold.

```bash
python run_test.py backbone
python run_test.py backbone "Express KRAS in HEK293 using pcDNA3.1 with Gibson"
```

The backbone name extracted by the LLM (e.g. `pcDNA3.1`) is searched first; host-derived defaults (pcDNA3.1 for mammalian, pUC19 for bacterial) are tried as fallback.

### Feature agent

Runs Intent then fetches all regulatory element sequences (promoter, polyA/terminator, selection marker, ori, tags) from NCBI and Addgene.

```bash
python run_test.py feature
python run_test.py feature "Express TP53 in HEK293 with CMV promoter and BGH polyA"
```

Feature resolution priority per element:
1. **Known fixed accession** — e.g. Ampicillin always fetches `V00613` (TEM-1, 800–900 bp), no search
2. **NCBI nuccore** — dynamic title/field search, annotated sub-feature extraction, length validation
3. **Addgene** — catalog search by `promoters=`, `resistance_marker=`, or `tags=` filter, GenBank download, annotation extraction

---

## Running the full pipeline

```bash
python -m plasmid_pipeline.run_pipeline
```

Or from Python:

```python
import asyncio
from plasmid_pipeline.orchestrator import PipelineOrchestrator

async def main():
    orch = PipelineOrchestrator()
    result = await orch.run("Express TP53 in HEK293 cells with CMV promoter, Gibson assembly")
    print(result.export_output.summary_json)

asyncio.run(main())
```

---

## Logs

All agent activity is logged to `agent_conversation.log` (or `$AGENT_LOG_PATH`).

Each line is timestamped and prefixed with the agent name:

```
2026-03-15 12:00:01 [INTENT] INPUT Express TP53 in HEK293...
2026-03-15 12:00:03 [INTENT] OUTPUT {"gene_symbol": "TP53", ...}
2026-03-15 12:00:03 [GENE] INPUT {"gene_symbol": "TP53", ...}
2026-03-15 12:00:05 [BACKBONE] INPUT {"expression_host": "HEK293", ...}
2026-03-15 12:00:05 [FEATURE] INPUT {"promoter": "CMV", ...}
...
```

Gene, Backbone, and Feature log entries will interleave because they run concurrently.

---

## Project structure

```
plasmid_pipeline/
├── models.py                   # all Pydantic input/output models
├── orchestrator.py             # pipeline coordinator (concurrent execution)
├── intent_agent.py             # LLM-based request parser (OpenAI)
├── gene_agent.py               # CDS fetcher (Ensembl + NCBI)
├── backbone_agent.py           # vector scaffold finder (NCBI)
├── feature_agent.py            # regulatory element fetcher (NCBI + Addgene)
├── plasmid_feature_extractor.py # Biopython GenBank annotation extractor
├── addgene_client.py           # Addgene API client
├── expression_agent.py         # host/codon compatibility checker
├── construct_agent.py          # cassette assembler
├── assembly_agent.py           # fragment + primer designer
├── export_agent.py             # GenBank / JSON / map exporter
├── web_utils.py                # async HTTP helpers (httpx)
└── api_server.py               # Flask REST API wrapper

run_test.py                     # per-agent test runner (intent → agent chain)
logging_utils.py                # shared file logger
requirements.txt
```

---

## Key models (`plasmid_pipeline/models.py`)

| Model | Used by |
|---|---|
| `IntentInput` / `IntentOutput` | IntentAgent |
| `GeneInput` / `GeneOutput` | GeneAgent |
| `BackboneInput` / `BackboneOutput` | BackboneAgent |
| `FeatureInput` / `FeatureOutput` | FeatureAgent |
| `ResolvedFeature` | FeatureAgent output; consumed by ConstructAgent |
| `ExpressionInput` / `ExpressionOutput` | ExpressionAgent |
| `ConstructInput` / `ConstructOutput` | ConstructAgent |
| `AssemblyInput` / `AssemblyOutput` | AssemblyAgent |
| `ExportInput` / `ExportOutput` | ExportAgent |
| `PipelineResult` | final orchestrator output |

`ResolvedFeature` fields: `name`, `type`, `sequence`, `source` (e.g. `"NCBI:V00613"` or `"Addgene:52535"`), `length`, `validated`, `position_hint`.
