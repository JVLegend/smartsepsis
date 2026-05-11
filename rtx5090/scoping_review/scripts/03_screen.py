#!/usr/bin/env python3
"""
Stage 3 — Title + abstract screening with two-pass LLM (Cohen kappa
agreement) and human-adjudication queue.

Backends (auto-detected, in order):
  1. Anthropic API     (env ANTHROPIC_API_KEY)
  2. Ollama local      (env OLLAMA_BASE_URL, default http://localhost:11434)

Two passes per record:
  pass_a — terse strict prompt
  pass_b — verbose contextual prompt
Both must independently judge include / exclude / unclear with reason.

Outputs:
  outputs/screening/decisions.jsonl   (per-record per-pass decisions)
  outputs/screening/included.jsonl    (both passes include)
  outputs/screening/excluded.jsonl    (both passes exclude)
  outputs/screening/unclear.jsonl     (disagreement OR any "unclear" -> human review queue)
  outputs/audit/screening_prompts.txt (exact prompts used)
  outputs/audit/screening_kappa.json  (Cohen's kappa between passes)

Checkpoint every 50 records. Resumable.
"""
from __future__ import annotations
import os, json, time, argparse, hashlib
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[1]
INP  = ROOT / "outputs" / "dedup" / "all_records.jsonl"
OUT  = ROOT / "outputs" / "screening"; OUT.mkdir(parents=True, exist_ok=True)
AUD  = ROOT / "outputs" / "audit"; AUD.mkdir(parents=True, exist_ok=True)
DONE = ROOT / "stage_03.done"

REVIEW_QUESTION = (
    "We are conducting a PRISMA-ScR scoping review on AI-augmented "
    "molecular diagnostics for ophthalmic infection. INCLUDE studies that "
    "report CRISPR-Cas-based or Nanopore metagenomic diagnostics applied "
    "to OCULAR samples or OCULAR pathogens (endophthalmitis, microbial "
    "keratitis, conjunctivitis, ocular surface infection). INCLUDE method "
    "papers if they explicitly position the assay for ocular samples. "
    "EXCLUDE: non-ocular contexts, non-CRISPR / non-Nanopore methods, "
    "narrative reviews without primary data, animal models without "
    "translational claims, and pre-2017 work (CRISPR-Dx emerged ~2017). "
    "INCLUDE if uncertain — mark as 'unclear' so a human can adjudicate."
)

PASS_A_PROMPT = (
    f"{REVIEW_QUESTION}\n\n"
    "Title: {title}\nAbstract: {abstract}\n\n"
    "Reply with strict JSON only: "
    '{{"decision":"include|exclude|unclear","reason":"one short sentence"}}'
)
PASS_B_PROMPT = (
    "You are a clinical epidemiologist screening abstracts for a scoping "
    f"review. {REVIEW_QUESTION}\n\n"
    "Title: {title}\nAbstract: {abstract}\n\n"
    "Think step by step about (i) population/sample, (ii) intervention/"
    "method, (iii) outcome reported, then decide. Output JSON only: "
    '{{"decision":"include|exclude|unclear","reason":"why"}}'
)

# Save prompts for the audit trail
(AUD / "screening_prompts.txt").write_text(
    "PASS A:\n" + PASS_A_PROMPT + "\n\n\nPASS B:\n" + PASS_B_PROMPT
)

def hash_id(rec: dict) -> str:
    h = hashlib.sha1()
    h.update((rec.get("title","") + rec.get("doi","") + (rec.get("pmid") or "")).encode("utf-8","ignore"))
    return h.hexdigest()[:16]

# -------- Backend abstraction --------------------------------------------
class LLM:
    def __init__(self):
        self.kind = None
        if os.environ.get("ANTHROPIC_API_KEY"):
            import anthropic; self.client = anthropic.Anthropic()
            self.kind = "anthropic"
        else:
            try:
                import ollama; self.client = ollama
                self.kind = "ollama"
            except ImportError:
                raise RuntimeError("No LLM backend available. Set ANTHROPIC_API_KEY or install ollama.")
        print(f"LLM backend: {self.kind}")

    def call(self, prompt: str, max_tokens: int = 200) -> str:
        if self.kind == "anthropic":
            msg = self.client.messages.create(
                model="claude-sonnet-4-6",
                max_tokens=max_tokens,
                messages=[{"role": "user", "content": prompt}],
            )
            return msg.content[0].text
        else:
            r = self.client.generate(model=os.environ.get("OLLAMA_MODEL","llama3.1:70b"), prompt=prompt, options={"temperature":0.1, "num_predict": max_tokens})
            return r["response"]

def parse_decision(raw: str) -> dict:
    # Find the first {...} JSON object
    import re
    m = re.search(r"\{[\s\S]*?\}", raw)
    if not m: return {"decision":"unclear","reason":f"unparseable: {raw[:120]}"}
    try:
        d = json.loads(m.group(0))
        d["decision"] = str(d.get("decision","unclear")).lower().strip()
        if d["decision"] not in ("include","exclude","unclear"):
            d["decision"] = "unclear"
        return d
    except Exception as e:
        return {"decision":"unclear","reason":f"json error: {e}"}

# -------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(); p.add_argument("--force", action="store_true")
    p.add_argument("--limit", type=int, default=0, help="limit number of records (for testing)")
    a = p.parse_args()
    if DONE.exists() and not a.force:
        print("Stage 3 already done."); return

    llm = LLM()
    decisions_path = OUT / "decisions.jsonl"
    # Resume: skip ids already decided
    seen = set()
    if decisions_path.exists():
        for ln in decisions_path.open():
            try: seen.add(json.loads(ln)["id"])
            except: pass
    print(f"Resume: {len(seen)} records already decided.")

    records = []
    for ln in INP.open():
        r = json.loads(ln); r["_id"] = hash_id(r)
        if r["_id"] not in seen: records.append(r)
    if a.limit: records = records[:a.limit]
    print(f"To screen: {len(records)}")

    out_dec = decisions_path.open("a")
    n_a_include = n_b_include = n_agree = 0
    for i, r in enumerate(records, 1):
        prompt_a = PASS_A_PROMPT.format(title=r.get("title",""), abstract=(r.get("abstract","") or "")[:4000])
        prompt_b = PASS_B_PROMPT.format(title=r.get("title",""), abstract=(r.get("abstract","") or "")[:4000])
        try:
            dec_a = parse_decision(llm.call(prompt_a))
            dec_b = parse_decision(llm.call(prompt_b))
        except Exception as e:
            print(f"  LLM error on {r['_id']}: {e}"); continue
        rec = {"id": r["_id"], "title": r.get("title",""), "doi": r.get("doi",""),
               "pmid": r.get("pmid"), "source": r.get("source",""),
               "pass_a": dec_a, "pass_b": dec_b}
        out_dec.write(json.dumps(rec) + "\n"); out_dec.flush()
        n_a_include += dec_a["decision"]=="include"
        n_b_include += dec_b["decision"]=="include"
        n_agree     += dec_a["decision"]==dec_b["decision"]
        if i % 50 == 0:
            print(f"  {i}/{len(records)} A_inc={n_a_include} B_inc={n_b_include} agree%={n_agree/i:.1%}")
    out_dec.close()

    # Reduce
    included = []; excluded = []; unclear = []
    for ln in decisions_path.open():
        d = json.loads(ln)
        a_dec, b_dec = d["pass_a"]["decision"], d["pass_b"]["decision"]
        if a_dec == "include" and b_dec == "include":
            included.append(d)
        elif a_dec == "exclude" and b_dec == "exclude":
            excluded.append(d)
        else:
            unclear.append(d)
    for name, arr in [("included", included), ("excluded", excluded), ("unclear", unclear)]:
        path = OUT / f"{name}.jsonl"
        with path.open("w") as f:
            for r in arr: f.write(json.dumps(r) + "\n")
        print(f"  {name}: {len(arr)} -> {path.relative_to(ROOT)}")

    # Cohen's kappa
    try:
        from sklearn.metrics import cohen_kappa_score
        a_arr = []; b_arr = []
        for ln in decisions_path.open():
            d = json.loads(ln)
            a_arr.append(d["pass_a"]["decision"]); b_arr.append(d["pass_b"]["decision"])
        kappa = cohen_kappa_score(a_arr, b_arr)
        (AUD / "screening_kappa.json").write_text(json.dumps({"cohen_kappa": kappa, "n": len(a_arr)}, indent=2))
        print(f"  Cohen's kappa (pass A vs pass B): {kappa:.3f}")
    except ImportError:
        pass

    DONE.write_text(time.strftime("%FT%TZ"))

if __name__ == "__main__":
    main()
