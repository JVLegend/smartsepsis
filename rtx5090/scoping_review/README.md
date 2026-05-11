# Scoping Review — Claude Code Orchestration Pack

**Audience:** Claude Code agent running on the user's RTX 5090 Linux box.
Working directory: this folder.

**Goal:** execute the PRISMA-ScR scoping review on **AI-augmented molecular
diagnostics for ophthalmic infection — CRISPR-Cas, Nanopore metagenomics,
and the open methodologic gaps**, producing a defensible included-studies
set + extracted data + draft figures and tables.

The full protocol, search strategy, and manuscript skeleton are already
written:

- `../../preprint/scoping_review_protocol.md`
- `../../preprint/scoping_review_search_strategy.md`
- `../../preprint/scoping_review_skeleton.md`

Read those first. This pack automates the mechanical execution.

---

## Execution plan (5 stages, each script idempotent / resumable)

| # | Stage | Script | Time | Output |
|---|---|---|---|---|
| 1 | Fetch raw search results from free APIs (PubMed, bioRxiv, medRxiv, ClinicalTrials.gov) | `scripts/01_fetch_searches.py` | 10–30 min | `outputs/raw/<source>.jsonl` |
| 2 | Deduplicate by DOI / fuzzy title | `scripts/02_dedup.py` | 1–2 min | `outputs/dedup/all_records.jsonl` |
| 3 | Title + abstract screening (LLM-assisted, two-reviewer simulation: Claude API or local Llama 3.1 70B) | `scripts/03_screen.py` | 30–60 min | `outputs/screening/{included,excluded,unclear}.jsonl` |
| 4 | Structured data extraction from included full-texts + abstracts | `scripts/04_extract.py` | 30–90 min | `outputs/extracted/charting_table.tsv` |
| 5 | Synthesis — PRISMA flow diagram numbers, tables 1–4, evidence-map JSON | `scripts/05_synthesize.py` | 2–5 min | `outputs/synthesis/` |

After stage 5, the user reviews `outputs/synthesis/` and drops the numbers
into `preprint/scoping_review_skeleton.md` to produce the final
manuscript.

---

## Two-reviewer constraint

PRISMA-ScR requires **two independent reviewers** for title/abstract
screening with a third for conflict resolution. We simulate this with
**two LLM passes using different models or different prompt phrasings**,
then any disagreement is flagged for human (you the user) review. This is
documented honestly in the methods — it is LLM-assisted screening with
human adjudication, not LLM-replacing-human screening.

---

## Pre-flight requirements

- Python 3.11 (re-use `../requirements.txt` venv)
- Internet access for E-utilities, bioRxiv API, ClinicalTrials.gov
- Either:
  - `ANTHROPIC_API_KEY` env var (preferred — uses Claude for screening); OR
  - Local model server (Ollama with `llama3.1:70b` or `qwen2.5:72b` on
    GPU). Set `OLLAMA_BASE_URL=http://localhost:11434` in env.
- No paid database access required for stage 1; we are upfront in the
  manuscript about which databases were and were not searched.

Embase / Scopus / Web of Science require institutional library access and
manual export. The protocol already accommodates this: for those bases,
the human exports `.ris` or `.csv` from the institutional portal and
drops files into `outputs/raw/manual/`. Stage 2 picks them up.

---

## How to run

```bash
cd rtx5090/scoping_review
source ../.venv/bin/activate      # reuse venv from rtx5090/ pack
pip install -r requirements.txt   # adds tldextract, rapidfuzz, etc.

export NCBI_EMAIL="iaparamedicos@gmail.com"
export ANTHROPIC_API_KEY="..."     # if available

# Stages
python scripts/01_fetch_searches.py
python scripts/02_dedup.py
python scripts/03_screen.py        # checkpoints every 50 records
python scripts/04_extract.py
python scripts/05_synthesize.py
```

Each stage writes a `stage_NN.done` sentinel; re-running skips completed
stages. Force re-run with `--force`.

---

## Honesty rules (carry into manuscript Methods)

The manuscript Methods section MUST disclose:
- LLM-assisted screening with [Claude 3.5 Sonnet OR Llama 3.1 70B];
  inter-pass agreement Cohen κ reported; human adjudication for all
  disagreements and a 10% random sample of agreements.
- Which databases were searched directly (PubMed, bioRxiv, medRxiv,
  ClinicalTrials.gov, WHO ICTRP) vs imported via institutional access
  (Embase, Scopus, WoS) vs not searched.
- Reproducibility: every stage's deterministic seed and the LLM prompts
  used are committed under `outputs/audit/`.

---

## Outputs to commit back to the repo

After stage 5, `99_pack_results.sh` (in the parent `rtx5090/`) will
include `scoping_review/outputs/synthesis/` in the tarball.

The user manually copies the final `charting_table.tsv` and synthesis
numbers into the manuscript skeleton.
