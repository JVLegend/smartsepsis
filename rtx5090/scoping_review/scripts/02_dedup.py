#!/usr/bin/env python3
"""
Stage 2 — Deduplicate by DOI (exact) then by title (rapidfuzz, ratio >= 92).
Reads all .jsonl under outputs/raw/ (including outputs/raw/manual/ for
institutional exports like Embase/Scopus/WoS .ris-converted .jsonl).
"""
from __future__ import annotations
import json, time, argparse
from pathlib import Path
from rapidfuzz import fuzz, process

ROOT = Path(__file__).resolve().parents[1]
RAW  = ROOT / "outputs" / "raw"
OUT  = ROOT / "outputs" / "dedup"
OUT.mkdir(parents=True, exist_ok=True)
DONE = ROOT / "stage_02.done"

def normalize_record(r: dict) -> dict:
    out = {"source": r.get("source",""), "raw": r}
    out["pmid"] = r.get("PMID") or r.get("pmid")
    out["doi"]  = (r.get("AID") or r.get("doi") or "").lower() if isinstance(r.get("AID"), str) else r.get("doi","")
    out["title"] = (r.get("TI") or r.get("title") or r.get("BriefTitle") or "").strip()
    out["abstract"] = r.get("AB") or r.get("abstract") or r.get("BriefSummary") or ""
    out["year"] = r.get("DP") or r.get("date") or r.get("StartDate") or ""
    out["url"] = (
        f"https://pubmed.ncbi.nlm.nih.gov/{out['pmid']}/" if out["pmid"]
        else f"https://www.biorxiv.org/content/{r.get('doi')}" if r.get("source") in ("biorxiv","medrxiv") and r.get("doi")
        else f"https://clinicaltrials.gov/study/{r.get('NCTId')}" if r.get("source") == "clinicaltrials" else ""
    )
    return out

def main():
    p = argparse.ArgumentParser(); p.add_argument("--force", action="store_true")
    a = p.parse_args()
    if (ROOT / "stage_02.done").exists() and not a.force:
        print("Stage 2 already done."); return

    records = []
    for src in sorted(RAW.rglob("*.jsonl")):
        with src.open() as f:
            for line in f:
                if line.strip(): records.append(normalize_record(json.loads(line)))
    print(f"Read {len(records)} raw records.")

    # DOI dedup
    seen_doi = {}
    by_doi = []
    for r in records:
        k = (r["doi"] or "").strip().lower()
        if k and k in seen_doi: continue
        if k: seen_doi[k] = True
        by_doi.append(r)
    print(f"After DOI dedup: {len(by_doi)}")

    # Title dedup with rapidfuzz
    kept = []; titles = []
    for r in by_doi:
        t = r["title"].lower().strip()
        if not t:
            kept.append(r); continue
        if titles:
            best = process.extractOne(t, titles, scorer=fuzz.token_set_ratio)
            if best and best[1] >= 92:
                continue
        kept.append(r); titles.append(t)
    print(f"After title dedup: {len(kept)}")

    out = OUT / "all_records.jsonl"
    with out.open("w") as f:
        for r in kept: f.write(json.dumps(r) + "\n")
    print(f"Wrote {out.relative_to(ROOT)}")

    (ROOT / "stage_02.done").write_text(time.strftime("%FT%TZ"))

if __name__ == "__main__":
    main()
