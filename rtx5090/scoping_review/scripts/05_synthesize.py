#!/usr/bin/env python3
"""
Stage 5 — Synthesis. Produces:
  - prisma_flow_numbers.json    (for the PRISMA-ScR Figure 1)
  - table1_included_studies.tsv
  - table2_target_conditions.tsv
  - table3_technology_stack.tsv
  - table4_validation_stage.tsv
  - evidence_map.json           (for the Figure 2 evidence map)
  - synthesis_notes.md          (human-readable summary)

Everything reads from outputs/screening/ and outputs/extracted/.
"""
from __future__ import annotations
import json, time, argparse, csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SC   = ROOT / "outputs" / "screening"
EX   = ROOT / "outputs" / "extracted"
OUT  = ROOT / "outputs" / "synthesis"; OUT.mkdir(parents=True, exist_ok=True)
RAW  = ROOT / "outputs" / "raw"
DONE = ROOT / "stage_05.done"

def count_jsonl(p):
    if not p.exists(): return 0
    return sum(1 for _ in p.open() if _.strip())

def main():
    p = argparse.ArgumentParser(); p.add_argument("--force", action="store_true"); a = p.parse_args()
    if DONE.exists() and not a.force: print("Stage 5 already done."); return

    # PRISMA numbers
    raw_total = sum(count_jsonl(f) for f in RAW.rglob("*.jsonl"))
    dedup_total = count_jsonl(ROOT/"outputs"/"dedup"/"all_records.jsonl")
    included = count_jsonl(SC/"included.jsonl")
    excluded = count_jsonl(SC/"excluded.jsonl")
    unclear  = count_jsonl(SC/"unclear.jsonl")
    flow = {
        "records_identified_from_databases": raw_total,
        "records_after_dedup": dedup_total,
        "records_screened": dedup_total,
        "records_included_after_screening_llm_pass": included,
        "records_excluded_after_screening_llm_pass": excluded,
        "records_unclear_for_human_adjudication": unclear,
    }
    (OUT/"prisma_flow_numbers.json").write_text(json.dumps(flow, indent=2))
    print(f"PRISMA numbers: {flow}")

    # Tables from extracted/charting_table.tsv
    charting = EX/"charting_table.tsv"
    if not charting.exists():
        print("No charting table yet — stage 4 incomplete."); return
    rows = list(csv.DictReader(charting.open(), delimiter="\t"))
    print(f"Extracted rows: {len(rows)}")

    # Table 1 — included studies (id, title, year, country, design, sample_type)
    with (OUT/"table1_included_studies.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id","title","year","country","design","sample_type"])
        for r in rows:
            w.writerow([r["id"], r["title"], r["year"], r["country"], r["design"], r["sample_type"]])

    # Table 2 — target conditions and pathogens
    pathogen_counts = Counter()
    for r in rows:
        for sp in (r["target_pathogens"] or "").split(","):
            sp = sp.strip()
            if sp: pathogen_counts[sp] += 1
    with (OUT/"table2_target_conditions.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["pathogen","studies_n"])
        for sp, n in pathogen_counts.most_common(): w.writerow([sp, n])

    # Table 3 — technology stack
    tech = Counter()
    for r in rows:
        tech[(r["cas_enzyme"] or "none", r["amplification_method"] or "none")] += 1
    with (OUT/"table3_technology_stack.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["cas_enzyme","amplification_method","studies_n"])
        for (c,a), n in tech.most_common(): w.writerow([c, a, n])

    # Table 4 — validation stage
    val = Counter(r["validation_stage"] or "unknown" for r in rows)
    with (OUT/"table4_validation_stage.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["validation_stage","studies_n"])
        for v, n in val.most_common(): w.writerow([v, n])

    # Evidence map (Figure 2 input): pathogen × validation_stage
    evidence = defaultdict(lambda: Counter())
    for r in rows:
        for sp in (r["target_pathogens"] or "").split(","):
            sp = sp.strip()
            if sp:
                evidence[sp][r["validation_stage"] or "unknown"] += 1
    (OUT/"evidence_map.json").write_text(json.dumps({k: dict(v) for k,v in evidence.items()}, indent=2))

    # Human-readable summary
    md = [
        "# Scoping Review — Synthesis Notes",
        f"_Generated: {time.strftime('%FT%TZ')}_",
        "",
        f"- Records identified: **{raw_total}**",
        f"- After dedup: **{dedup_total}**",
        f"- Included after LLM screening (both passes): **{included}**",
        f"- Excluded: **{excluded}**",
        f"- Unclear (human adjudication queue): **{unclear}**",
        "",
        "## Pathogens covered (top 10)",
    ]
    for sp, n in pathogen_counts.most_common(10): md.append(f"- {sp} — {n} studies")
    md.append("\n## Technology stack")
    for (c,a), n in tech.most_common(): md.append(f"- {c} + {a} — {n}")
    md.append("\n## Validation stage")
    for v, n in val.most_common(): md.append(f"- {v} — {n}")
    md.append("\n## Open gap (analyst note)")
    md.append("Combine the evidence map (`evidence_map.json`) with Table 4 to identify (pathogen × stage) cells that are empty. These are the scientific gaps the SmartSepsis-Oph program targets.")
    (OUT/"synthesis_notes.md").write_text("\n".join(md))
    print(f"Wrote synthesis to {OUT.relative_to(ROOT)}")

    DONE.write_text(time.strftime("%FT%TZ"))

if __name__ == "__main__":
    main()
