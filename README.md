# Plasmid Creation Agent

A natural-language plasmid design tool. Describe the construct you want in plain English and the pipeline interprets your intent, assembles the sequence from curated datasets (with automatic web fallback), and exports a ready-to-open GenBank file.

---

## Why this exists

Designing a plasmid still requires opening multiple browser tabs — looking up promoter sequences on Addgene, copying RBS variants from the iGEM registry, checking polyA signals in the literature — then manually stitching everything together in SnapGene or Benchling. This is slow and error-prone.

This app removes that friction. You type one sentence (`"Express EGFP in HEK293 cells with CMV promoter"`), and the pipeline:

1. Parses your intent and fills in all missing defaults automatically
2. Looks up every element (promoter, RBS, linker, polyA, MCS, backbone…) from a local curated database
3. Falls back to NCBI / Addgene for anything not in the database — fetching clean sequences, not whole plasmids
4. Assembles the full construct and exports a `.gbk` file you can open directly in SnapGene or Benchling

---

## Setup

### 1. Clone & enter the project

```bash
cd plasmid-creation-agent-2
```

### 2. Create and activate the virtual environment

```bash
# Create (only needed once)
python -m venv venv

# Activate — Windows
venv\Scripts\activate

# Activate — macOS / Linux
source venv/bin/activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Configure API keys

Create a `.env` file in the project root (or edit the existing one):

```
OPENAI_API_KEY=sk-...        # required — used by the Intent Agent
ADDGENE_API_KEY=devapi-...   # optional — improves backbone retrieval
```

---

## Running the app

### Frontend (browser UI)

```bash
python app_server.py
```

Open **http://localhost:5000** in your browser.

### CLI / quick tests

```bash
# Parse intent only (fast, no web calls)
python run_test.py intent "Express EGFP in HEK293 cells"

# Intent + construction (resolves all elements)
python run_test.py construction "Express EGFP in HEK293 cells"

# Full pipeline + GBK preview
python run_test.py pipeline "Express EGFP in HEK293 cells"

# Full pipeline + GBK file written to disk
python run_test.py output "Express EGFP in HEK293 cells"
```

Default request (used when no argument is given):
`"Create a fusion plasmid with Gene A TP53 and Gene B EGFR in HEK293 using CMV promoter"`

---

## The three agents

### 1. Intent Agent (`pipeline_agent/intent_agent.py`)

**What it does:** Converts a free-text user request into a structured plasmid specification.

Uses the OpenAI API (GPT-4o) to extract:

| Field | Example |
|---|---|
| `gene_symbol` | `EGFP` |
| `expression_host` | `HEK293` |
| `backbone` | `pcDNA3.1` |
| `promoter` | `CMV` |
| `rbs` | `Kozak Consensus` (or `RBS_B0034` for E. coli) |
| `linker_sequence` | `GGGGSx3` |
| `terminator` / `polyA` | `BGH polyA` |
| `selection_marker` | `Neomycin (G418)` |
| `construct_type` | `standard` \| `fusion` \| `multi_cassette_fusion` |
| `cloning_site` | `pUC19 MCS` |
| … and more | |

Automatically fills in sensible defaults — if the user doesn't specify a promoter for a mammalian construct, it defaults to CMV; for E. coli it defaults to T7.

---

### 2. Construction Agent (`pipeline_agent/plasmid_construction_agent.py`)

**What it does:** Builds the full plasmid sequence element by element in the correct biological order.

Resolution strategy (in priority order):

```
1. Local dataset (plasmid_*.csv files)  →  instant, no network
2. Direct NCBI accession fetch           →  when CSV row has a known accession
3. Addgene API (by name or ID)           →  full plasmid sequences
4. NCBI text search (size-filtered)      →  element-type-aware, avoids whole-plasmid returns
5. Placeholder                           →  last resort, flagged with a warning
```

Key intelligence built in:
- **Host detection** — mammalian vs. bacterial changes which slots are included (RBS omitted for mammalian; Kozak added instead)
- **Backbone feature detection** — scans the backbone sequence for promoters, polyA signals, resistance markers, and MCS; skips adding elements the backbone already has
- **MCS auto-selection** — picks the correct polylinker for the backbone (pET MCS for pET vectors, pcDNA3.1 MCS for mammalian, pUC19 MCS otherwise); skips entirely if backbone already has ≥3 restriction sites
- **Cargo excision** — if the fetched backbone contains a reporter gene (GFP, mCherry…), it is automatically removed before inserting the new gene of interest
- **Clean CDS fetch** — for cargo genes not in the database, uses NCBI `fasta_cds_nt` to return only the coding sequence, not an entire plasmid record

---

### 3. Output Agent (`pipeline_agent/output_agent.py`)

**What it does:** Converts the assembled construct into a GenBank (`.gbk`) file.

- Builds a `BioPython` `SeqRecord` from the assembled sequence
- Annotates every element as a named feature with its type and source
- Writes a standard GenBank flat file — directly openable in SnapGene, Benchling, ApE, or any sequence viewer
- Passes through any warnings from the construction step (e.g. "element fetched from web — verify sequence")

---

## Dataset (`data/`)

All local sequence data lives in `data/plasmid_*.csv`. Every file is loaded automatically at startup.

| File | Contents | Key columns |
|---|---|---|
| `plasmid_Backbones.csv` | 18 common vectors (pUC19, pET-28a, pcDNA3.1, pAAV, pLenti…) | NCBI Accession, Addgene ID |
| `plasmid_Promoters.csv` | 8 promoters — T7, araBAD, tac, CMV, CAG… | DNA Sequence, Strength, Host |
| `plasmid_RBS.csv` | 9 ribosome binding sites — B0034, Shine-Dalgarno, Kozak, T7 RBS… | DNA Sequence, Strength |
| `plasmid_Linkers.csv` | 10 linker variants — GGGGS×3, GGS, EAAAK, PG, helical… | DNA Sequence, Amino Acid Sequence |
| `plasmid_Terminators.csv` | 5 bacterial terminators — BBa_B0010, B0012, T7, rho-dep, HDV | DNA Sequence, Efficiency |
| `plasmid_PolyA.csv` | 9 polyA signals — BGH, SV40, hGH, SPA, TK… | DNA Sequence, NCBI Accession |
| `plasmid_MCS.csv` | 8 cloning sites — pUC19, pBluescript, pET, pcDNA3.1, Golden Gate… | DNA Sequence, Restriction Sites |
| `plasmid_Origins.csv` | 10 origins — ColE1, p15A, pSC101, f1, SV40 ori, 2μ… | NCBI Accession, Copy Number |
| `plasmid_Cargo.csv` | 10 reporters — EGFP (full seq), mCherry, Firefly Luc, NanoLuc… | DNA Sequence, Amino Acid Sequence |
| `plasmid_Tags.csv` | Affinity and epitope tags (His, FLAG, HA, Strep, MBP…) | DNA Sequence |
| `plasmid_NLS.csv` | Nuclear localisation signals (SV40 NLS, bipartite…) | DNA Sequence |
| `plasmid_Signal_Peptides.csv` | Signal peptides for secretion | DNA Sequence |
| `plasmid_Linkers.csv` | Protein fusion linkers | DNA Sequence |
| `plasmid_Fusion_Proteins.csv` | Common fusion partners (GST, MBP, SUMO…) | DNA Sequence |
| `plasmid_Protease_Cleavage_Sites.csv` | TEV, thrombin, Factor Xa, 3C… | DNA Sequence |

### How the database is used

When the construction agent needs an element (e.g. "BGH polyA"):

1. It scans **all** `plasmid_*.csv` files for a name match
2. If the row has a `DNA Sequence` column → used directly
3. If the row has an `Amino Acid Sequence` but no DNA → translated using a standard codon table
4. If the row has an `NCBI Accession` but no sequence → the accession is used for a direct NCBI fetch (guaranteed correct record, no guessing)
5. If nothing is found locally → web search with element-type-aware size filters (e.g. polyA signals searched in the 150–700 bp range to avoid returning full plasmids)

### Adding your own elements

Add a row to the appropriate CSV. The minimum required columns are:

```
Element Name   (or Linker Name / Terminator Name / Backbone Name)
DNA Sequence (5' → 3')
```

Everything else is optional metadata. The agent picks up changes automatically on the next run — no code changes needed.
