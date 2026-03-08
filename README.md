# Plasmid Pipeline Multi-Agent System (Python)

## Overview

**Plasmid Pipeline** is a multi-agent biological design system written in Python. It transforms *imprecise, natural-language* requests into precise, machine-readable build plans for automated plasmid construction.

---

## What is Plasmid Pipeline?

The pipeline acts as an intelligent layer between the user and the plasmid design logic. You can input complex biological requests, and the system automatically generates detailed plans suitable for laboratory automation or synthetic biology work.

### Example: Natural Language to Machine Plan

> “Build a TP53 expression plasmid for HEK cells with a CMV promoter using Gibson assembly.”

<p align="center">
  <img src="public/1.jpg" width="420" alt="Example UI Input"/>
</p>

This input is automatically converted into structured data for downstream tools:

<p align="center">
  <img src="public/3.jpg" width="420" alt="Example Output Conversion"/>
  <img src="public/4.jpg" width="420" alt="Example Output Conversion"/>
</p>

---

## Key Features

- **Requirement Extraction**: Identifies gene, host, promoter, and assembly method from free text.
- **Backbone Selection**: Suggests optimal vectors for your requirements.
- **Construct Briefing**: Generates plasmid maps and primer designs.
- **Data Export**: Produces JSON outputs ready for lab automation.

<p align="center">
  <img src="public/2.jpg" width="400" alt="Pipeline output overview"/>
</p>

---

# Getting Started

### 1. Python Version

> **Requires Python 3.10 or higher (3.11/3.12 recommended). Use a clean virtual environment.**

### 2. Virtual Environment Setup

#### macOS / Linux

```bash
cd a2a-python-test
python3 -m venv venv
source venv/bin/activate
# To deactivate: deactivate
```

#### Windows (PowerShell)

```powershell
cd a2a-python-test
python -m venv venv
.\venv\Scripts\Activate.ps1
# To deactivate: deactivate
```

> *If PowerShell blocks activation, run:*  
> `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`

### 3. Install Dependencies

```bash
pip install -r requirements.txt
```
*Installs*: `a2a-sdk`, `python-dotenv`, `openai`, `uvicorn`, `httpx`, etc.

### 4. Set Environment Variables

You **must** supply an OpenAI API key:

- **macOS / Linux**  
  `export OPENAI_API_KEY="sk-..."`
- **Windows (Powershell)**  
  `$env:OPENAI_API_KEY="sk-..."`
- **Or**: Create a `.env` file in the root directory:  
  `OPENAI_API_KEY=sk-...`

---

## Project Architecture

- **Orchestrator Agent**: Coordinates the workflow across agents.
- **Specialized Agents**: Handle individual tasks (design, backbone, QC, etc.).

### Agent Servers

> Open **four separate terminals**, then run each:

#### Terminal 1 — Orchestrator Agent

```bash
python -m agents.orchestrator_agent.main --port 9090
```

#### Terminal 2 — Design Agent

```bash
python -m agents.design_agent.main --port 9101
```

#### Terminal 3 — Backbone Agent

```bash
python -m agents.backbone_agent.main --port 9102
```

#### Terminal 4 — QC Agent

```bash
python -m agents.risk_agent.main --port 9103
```

> **TIP:**  
> Rename any legacy module names (like `crypto_agent`) to biological names (e.g., `design_agent`) for readability.

---

## Running the Client

Open a **fifth terminal** and send your request:

```bash
python -m client.main --card-url http://127.0.0.1:9090 --log agent_conversation.log
```

---

## Example JSON Output

```json
{
  "gene": "TP53",
  "host": "mammalian",
  "promoter": "CMV",
  "vector_constraint": "plasmid",
  "assembly_preference": "gibson",
  "output": "genbank"
}
```

---

## Logs and Agent Conversation

Check `agent_conversation.log` to follow the agent workflow:

- `client → orchestrator`
- `orchestrator → sub-agent`
- `sub-agent → orchestrator`
- `orchestrator → client`

---
