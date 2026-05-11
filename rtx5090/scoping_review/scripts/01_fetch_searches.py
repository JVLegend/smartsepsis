#!/usr/bin/env python3
"""
Stage 1 — Fetch raw search results from free APIs.

Sources:
  - PubMed (E-utilities, free, no key needed; respect 3 req/s)
  - bioRxiv (free JSON API)
  - medRxiv (free JSON API)
  - ClinicalTrials.gov (free v2 API)
  - WHO ICTRP (HTML/CSV, optional — disabled by default; manual export OK)

Each source writes outputs/raw/<source>.jsonl. Re-running skips sources
already complete unless --force.
"""
from __future__ import annotations
import os, sys, json, time, argparse, urllib.parse
from pathlib import Path
import requests
from Bio import Entrez

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "outputs" / "raw"
OUT.mkdir(parents=True, exist_ok=True)
DONE = ROOT / "stage_01.done"

# -----------------------------------------------------------------------
# QUERIES
# Derived from preprint/scoping_review_search_strategy.md. We keep these
# inline so this script is self-documenting and trivially editable.
# -----------------------------------------------------------------------
QUERIES = {
    "pubmed_main": (
        '('
        '"CRISPR-Cas Systems"[MeSH] OR CRISPR[tiab] OR Cas12[tiab] OR Cas13[tiab] '
        'OR SHERLOCK[tiab] OR DETECTR[tiab] OR "collateral cleavage"[tiab]'
        ') AND ('
        'Nanopore[tiab] OR "Nanopore Sequencing"[MeSH] OR metagenomic*[tiab]'
        ' OR "Diagnostic Tests, Routine"[MeSH] OR "point-of-care"[tiab]'
        ' OR isothermal[tiab] OR RPA[tiab] OR LAMP[tiab]'
        ') AND ('
        'eye[tiab] OR ocular[tiab] OR ophthalm*[tiab] OR endophthalmitis[MeSH]'
        ' OR endophthalmitis[tiab] OR keratitis[MeSH] OR keratitis[tiab]'
        ' OR uveitis[tiab] OR conjunctivitis[tiab] OR vitreous[tiab]'
        ' OR cornea[tiab] OR "aqueous humor"[tiab]'
        ')'
    ),
    "pubmed_adjacent_nanopore": (
        '(Nanopore[tiab] OR "Nanopore Sequencing"[MeSH] OR metagenomic*[tiab])'
        ' AND (eye[tiab] OR ocular[tiab] OR ophthalm*[tiab]'
        ' OR endophthalmitis[tiab] OR keratitis[tiab] OR vitreous[tiab])'
    ),
    "pubmed_adjacent_crispr_dx": (
        '(CRISPR[tiab] OR Cas12[tiab] OR Cas13[tiab] OR SHERLOCK[tiab] OR DETECTR[tiab])'
        ' AND (diagnostic*[tiab] OR detection[tiab] OR "point-of-care"[tiab])'
        ' AND (bacterial[tiab] OR fungal[tiab] OR "infectious disease"[tiab])'
    ),
    "biorxiv": "CRISPR ophthalmic OR Nanopore endophthalmitis OR Cas12 ocular",
    "medrxiv": "CRISPR ophthalmic OR Nanopore endophthalmitis OR Cas12 ocular",
    "clinicaltrials": "CRISPR ophthalmic OR Nanopore endophthalmitis",
}

# -----------------------------------------------------------------------
def jsonl_write(path: Path, records):
    with path.open("w", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"  wrote {len(records):>5} records -> {path.relative_to(ROOT)}")

# -----------------------------------------------------------------------
# PubMed (E-utilities). Free, requires email; respect rate limits.
# -----------------------------------------------------------------------
def fetch_pubmed(name: str, query: str, retmax: int = 1000):
    Entrez.email = os.environ.get("NCBI_EMAIL", "anonymous@example.org")
    out = OUT / f"pubmed_{name}.jsonl"
    if out.exists() and not FORCE:
        print(f"  skip {out.name} (exists)"); return
    print(f"[PubMed] '{name}' …")
    h = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    ids = Entrez.read(h)["IdList"]
    print(f"  {len(ids)} PMIDs")
    records = []
    for i in range(0, len(ids), 200):
        chunk = ids[i:i+200]
        h = Entrez.efetch(db="pubmed", id=",".join(chunk), rettype="medline", retmode="text")
        text = h.read()
        # Trivial MEDLINE parser → records
        for blk in text.split("\n\n"):
            if not blk.strip(): continue
            rec = {"source": "pubmed", "query_name": name}
            cur_key = None
            for line in blk.split("\n"):
                if line.startswith("      "):  # continuation
                    if cur_key: rec[cur_key] = rec.get(cur_key, "") + " " + line.strip()
                else:
                    parts = line.split("- ", 1)
                    if len(parts) == 2:
                        k, v = parts[0].strip(), parts[1].strip()
                        cur_key = k
                        # multi-valued tags
                        if k in ("AU","AD","MH","FAU","OT","PT"):
                            rec.setdefault(k, []).append(v)
                        else:
                            rec[k] = v
            if rec.get("PMID"): records.append(rec)
        time.sleep(0.34)  # NCBI: 3 req/sec without key
    jsonl_write(out, records)

# -----------------------------------------------------------------------
# bioRxiv / medRxiv via JSON API
#   https://api.biorxiv.org/details/{server}/{from}/{to}/{cursor}
#   We use the public details endpoint and filter client-side because the
#   `q=` parameter is not officially supported.
# -----------------------------------------------------------------------
def fetch_rxiv(server: str, query: str, start_yyyy_mm_dd="2017-01-01"):
    out = OUT / f"{server}.jsonl"
    if out.exists() and not FORCE: print(f"  skip {out.name}"); return
    print(f"[{server}] crawling from {start_yyyy_mm_dd} …")
    terms = [t.strip().lower() for t in query.replace(" OR ", "|").split("|")]
    today = time.strftime("%Y-%m-%d")
    cursor = 0; records = []
    while True:
        url = f"https://api.biorxiv.org/details/{server}/{start_yyyy_mm_dd}/{today}/{cursor}"
        r = requests.get(url, timeout=30)
        if r.status_code != 200: break
        d = r.json()
        coll = d.get("collection", [])
        if not coll: break
        for it in coll:
            text = (it.get("title","") + " " + it.get("abstract","")).lower()
            if any(t.replace(" ", "") in text.replace(" ", "") for t in terms):
                records.append({"source": server, **it})
        cursor += len(coll)
        if cursor > 50000: break  # safety cap
        time.sleep(0.5)
    jsonl_write(out, records)

# -----------------------------------------------------------------------
# ClinicalTrials.gov v2
# -----------------------------------------------------------------------
def fetch_clinicaltrials(query: str):
    out = OUT / "clinicaltrials.jsonl"
    if out.exists() and not FORCE: print(f"  skip {out.name}"); return
    print("[ClinicalTrials.gov] …")
    base = "https://clinicaltrials.gov/api/v2/studies"
    params = {"query.term": query, "pageSize": 100, "fields": "NCTId,BriefTitle,OfficialTitle,Phase,Condition,InterventionName,BriefSummary,OverallStatus,StartDate"}
    records = []
    nxt = None
    while True:
        if nxt: params["pageToken"] = nxt
        r = requests.get(base, params=params, timeout=30)
        if r.status_code != 200: print(f"  HTTP {r.status_code}"); break
        d = r.json()
        for s in d.get("studies", []):
            records.append({"source":"clinicaltrials", **s})
        nxt = d.get("nextPageToken")
        if not nxt: break
        time.sleep(0.5)
    jsonl_write(out, records)

# -----------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    global FORCE; FORCE = a.force

    if DONE.exists() and not FORCE:
        print("Stage 1 already done. Use --force to re-run."); return

    fetch_pubmed("main", QUERIES["pubmed_main"])
    fetch_pubmed("adjacent_nanopore", QUERIES["pubmed_adjacent_nanopore"], retmax=500)
    fetch_pubmed("adjacent_crispr_dx", QUERIES["pubmed_adjacent_crispr_dx"], retmax=500)
    fetch_rxiv("biorxiv", QUERIES["biorxiv"])
    fetch_rxiv("medrxiv", QUERIES["medrxiv"])
    fetch_clinicaltrials(QUERIES["clinicaltrials"])

    DONE.write_text(time.strftime("%FT%TZ"))
    print("\nStage 1 done.")

if __name__ == "__main__":
    main()
