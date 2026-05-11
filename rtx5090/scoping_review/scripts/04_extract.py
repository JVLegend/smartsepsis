#!/usr/bin/env python3
"""
Stage 4 — Structured data extraction from included records.

For each record in outputs/screening/included.jsonl, ask the LLM to
extract a structured JSON object matching the PRISMA-ScR charting
template:

  study_id, year, country, design, sample_type, cas_enzyme,
  amplification_method, target_pathogens, target_genes,
  lod_reported, sensitivity_pct, specificity_pct, time_to_result_min,
  validation_stage, key_findings, limitations

Optionally pulls full-text via Unpaywall (if email + free PDF) for
better extraction. Otherwise uses title + abstract only.

Outputs:
  outputs/extracted/charting_table.tsv
  outputs/extracted/raw_extractions.jsonl
"""
from __future__ import annotations
import os, json, time, argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INC  = ROOT / "outputs" / "screening" / "included.jsonl"
RAW  = ROOT / "outputs" / "dedup" / "all_records.jsonl"
OUT  = ROOT / "outputs" / "extracted"; OUT.mkdir(parents=True, exist_ok=True)
DONE = ROOT / "stage_04.done"

EXTRACTION_PROMPT = """You are charting data for a PRISMA-ScR scoping review.

Title: {title}
Abstract: {abstract}

Extract the following fields as JSON. Use null if not reported. Be precise; do not invent.

{{
  "year": <int or null>,
  "country": "<str or null>",
  "design": "<original article | methods paper | review | case series | trial protocol | other>",
  "sample_type": "<vitreous tap | corneal scrape | aqueous humor | ocular surface swab | cultured isolate | other>",
  "cas_enzyme": "<Cas12a | Cas13 | Cas9 | none>",
  "amplification_method": "<RPA | LAMP | PCR | none>",
  "target_pathogens": ["list","of","species"],
  "target_genes": ["list","of","genes"],
  "lod_reported_copies_per_ml": <float or null>,
  "sensitivity_pct": <float or null>,
  "specificity_pct": <float or null>,
  "time_to_result_min": <float or null>,
  "validation_stage": "<in silico | analytical | clinical pilot | full clinical | regulatory>",
  "key_findings": "<one sentence>",
  "limitations": "<one sentence>"
}}

Respond with the JSON object only."""

class LLM:
    def __init__(self):
        if os.environ.get("ANTHROPIC_API_KEY"):
            import anthropic; self.client = anthropic.Anthropic(); self.kind = "anthropic"
        else:
            import ollama; self.client = ollama; self.kind = "ollama"
        print(f"LLM backend: {self.kind}")
    def call(self, prompt: str) -> str:
        if self.kind == "anthropic":
            r = self.client.messages.create(model="claude-sonnet-4-6", max_tokens=800,
                messages=[{"role":"user","content":prompt}])
            return r.content[0].text
        return self.client.generate(model=os.environ.get("OLLAMA_MODEL","llama3.1:70b"), prompt=prompt, options={"temperature":0.0,"num_predict":800})["response"]

def parse_json(raw: str) -> dict:
    import re
    m = re.search(r"\{[\s\S]+\}", raw)
    if not m: return {}
    try: return json.loads(m.group(0))
    except: return {}

def main():
    p = argparse.ArgumentParser(); p.add_argument("--force", action="store_true"); a = p.parse_args()
    if DONE.exists() and not a.force: print("Stage 4 already done."); return

    # Build id -> record map from dedup file (for abstract content)
    id2rec = {}
    import hashlib
    for ln in RAW.open():
        r = json.loads(ln)
        h = hashlib.sha1(((r.get("title","")+r.get("doi","")+(r.get("pmid") or "")).encode("utf-8","ignore"))).hexdigest()[:16]
        id2rec[h] = r
    llm = LLM()
    extractions = []
    seen = set()
    raw_path = OUT / "raw_extractions.jsonl"
    if raw_path.exists():
        for ln in raw_path.open():
            try: seen.add(json.loads(ln)["id"])
            except: pass
    out = raw_path.open("a")
    for ln in INC.open():
        d = json.loads(ln)
        if d["id"] in seen: continue
        rec = id2rec.get(d["id"], d)
        prompt = EXTRACTION_PROMPT.format(title=rec.get("title","")[:500], abstract=(rec.get("abstract","") or "")[:6000])
        try:
            raw = llm.call(prompt)
            data = parse_json(raw)
        except Exception as e:
            data = {"_error": str(e)}
        record = {"id": d["id"], "title": rec.get("title",""), "doi": rec.get("doi",""), "extracted": data}
        out.write(json.dumps(record) + "\n"); out.flush()
        extractions.append(record)
    out.close()

    # Build TSV table
    import csv
    cols = ["id","title","year","country","design","sample_type","cas_enzyme","amplification_method","target_pathogens","target_genes","lod_reported_copies_per_ml","sensitivity_pct","specificity_pct","time_to_result_min","validation_stage","key_findings","limitations"]
    with (OUT / "charting_table.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for ln in raw_path.open():
            r = json.loads(ln); e = r.get("extracted", {})
            row = [r["id"], r["title"][:120]]
            for c in cols[2:]:
                v = e.get(c)
                if isinstance(v, list): v = ", ".join(map(str, v))
                row.append(v if v is not None else "")
            w.writerow(row)
    print(f"Wrote {OUT.relative_to(ROOT)}/charting_table.tsv")
    DONE.write_text(time.strftime("%FT%TZ"))

if __name__ == "__main__":
    main()
