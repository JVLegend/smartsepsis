#!/usr/bin/env python3
"""
Análise de conservação: verifica se os guides atuais cobrem as variantes.
Para cada gene family, compara o spacer do guide contra as sequências das variantes.
"""

import os
import csv
import time
import requests

from config import BASE_DIR, SEQUENCES_DIR, GUIDES_DIR, REPORTS_DIR, NCBI_BASE_URL, NCBI_EMAIL
from utils import parse_fasta, reverse_complement


def fetch_variant_sequence(accession: str, name: str) -> str | None:
    """Baixa sequência de uma variante do NCBI."""
    url = f"{NCBI_BASE_URL}/efetch.fcgi"
    params = {
        "db": "nucleotide", "id": accession,
        "rettype": "fasta", "retmode": "text", "email": NCBI_EMAIL,
    }
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        content = resp.text.strip()
        if content.startswith(">"):
            lines = content.split("\n")
            seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
            return seq.upper()
    except Exception as e:
        print(f"    Erro ao buscar {name}: {e}")
    return None


def find_spacer_in_sequence(spacer: str, sequence: str) -> dict:
    """Procura o spacer (e seu rev complement) na sequência da variante."""
    spacer_up = spacer.upper()
    seq_up = sequence.upper()
    rc_spacer = reverse_complement(spacer_up)

    # Match exato
    if spacer_up in seq_up:
        pos = seq_up.index(spacer_up)
        return {"match": "exact", "mismatches": 0, "position": pos, "strand": "+"}
    if rc_spacer in seq_up:
        pos = seq_up.index(rc_spacer)
        return {"match": "exact", "mismatches": 0, "position": pos, "strand": "-"}

    # Match com mismatches (sliding window)
    best_mm = len(spacer_up)
    best_pos = -1
    best_strand = "+"

    for strand, query in [("+", spacer_up), ("-", rc_spacer)]:
        for i in range(len(seq_up) - len(query) + 1):
            window = seq_up[i:i + len(query)]
            mm = sum(1 for a, b in zip(query, window) if a != b)
            if mm < best_mm:
                best_mm = mm
                best_pos = i
                best_strand = strand

    if best_mm <= 3:
        return {"match": "near", "mismatches": best_mm, "position": best_pos, "strand": best_strand}
    elif best_mm <= 5:
        return {"match": "partial", "mismatches": best_mm, "position": best_pos, "strand": best_strand}
    else:
        return {"match": "none", "mismatches": best_mm, "position": best_pos, "strand": best_strand}


def load_guide_spacer(gene_family: str) -> str | None:
    """Carrega o melhor spacer de um gene family."""
    # Mapear variantes para o gene family base
    family_map = {
        "mecA": "mecA", "mecA1": "mecA", "mecA2": "mecA",
        "blaKPC": "blaKPC", "blaNDM": "blaNDM", "vanA": "vanA",
        "mcr": "mcr-1", "mcr-1": "mcr-1", "mcr-5": "mcr-1",
        "blaOXA-48": "blaOXA-48", "blaOXA-181": "blaOXA-48", "blaOXA-232": "blaOXA-48",
        "blaVIM": "blaVIM", "blaIMP": "blaIMP", "blaGES": "blaGES",
        "blaCTX-M": "blaCTX-M-15", "qnrS": "qnrS", "armA": "armA",
    }
    base = family_map.get(gene_family, gene_family)

    guides_path = os.path.join(GUIDES_DIR, f"{base}_cas12a_guides.tsv")
    if not os.path.exists(guides_path):
        return None
    with open(guides_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            return row["spacer_seq"]
    return None


def run_conservation_analysis(variants_csv: str = "targets_brazil_variants.csv"):
    """Roda análise de conservação para todas as variantes."""
    print("=" * 70)
    print("ANÁLISE DE CONSERVAÇÃO - Cobertura dos guides atuais")
    print("=" * 70)

    os.makedirs(REPORTS_DIR, exist_ok=True)

    # Carregar variantes
    variants = []
    with open(os.path.join(BASE_DIR, variants_csv)) as f:
        reader = csv.DictReader(f)
        for row in reader:
            variants.append(row)

    print(f"\nVariantes a analisar: {len(variants)}")

    # Agrupar por gene family
    families = {}
    for v in variants:
        fam = v.get("gene_family", v["name"].split("-")[0] if "-" in v["name"] else v["name"])
        if fam not in families:
            families[fam] = []
        families[fam].append(v)

    results = []
    summary = {}

    for family, family_variants in families.items():
        spacer = load_guide_spacer(family)
        if not spacer:
            print(f"\n[{family}] Sem guide disponível, pulando")
            continue

        print(f"\n{'─'*50}")
        print(f"[{family}] Spacer: {spacer} ({len(family_variants)} variantes)")

        covered = 0
        near = 0
        partial = 0
        missed = 0

        for v in family_variants:
            name = v["name"]
            acc = v["gene_accession"]

            # Buscar sequência
            seq = fetch_variant_sequence(acc, name)
            if not seq:
                print(f"  ✗ {name}: falha no fetch")
                results.append({"variant": name, "family": family, "accession": acc,
                                "match": "fetch_failed", "mismatches": -1})
                continue

            time.sleep(0.4)

            # Procurar spacer
            match = find_spacer_in_sequence(spacer, seq)

            if match["match"] == "exact":
                covered += 1
                icon = "✓"
            elif match["match"] == "near":
                near += 1
                icon = "~"
            elif match["match"] == "partial":
                partial += 1
                icon = "?"
            else:
                missed += 1
                icon = "✗"

            print(f"  {icon} {name:20s} {match['match']:8s} ({match['mismatches']}mm) pos={match['position']} {match['strand']}")
            results.append({
                "variant": name, "family": family, "accession": acc,
                "match": match["match"], "mismatches": match["mismatches"],
                "position": match["position"], "strand": match["strand"],
            })

        total = covered + near + partial + missed
        coverage_exact = (covered / total * 100) if total > 0 else 0
        coverage_near = ((covered + near) / total * 100) if total > 0 else 0

        summary[family] = {
            "total": total, "exact": covered, "near": near,
            "partial": partial, "missed": missed,
            "coverage_exact": round(coverage_exact, 1),
            "coverage_near": round(coverage_near, 1),
        }

        print(f"  Cobertura: {covered}/{total} exato ({coverage_exact:.0f}%) | +{near} near ({coverage_near:.0f}% com ≤3mm)")

    # Salvar relatório
    report_path = os.path.join(REPORTS_DIR, "conservation_analysis.tsv")
    with open(report_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["variant", "family", "accession", "match", "mismatches", "position", "strand"], delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    # Resumo
    print(f"\n{'='*70}")
    print("RESUMO DE COBERTURA")
    print(f"{'='*70}")
    print(f"{'Family':<15} {'Total':<8} {'Exato':<8} {'≤3mm':<8} {'Cob.Exata':<12} {'Cob.≤3mm':<12}")
    print("-" * 65)
    for fam, s in sorted(summary.items()):
        print(f"{fam:<15} {s['total']:<8} {s['exact']:<8} {s['exact']+s['near']:<8} {s['coverage_exact']:.0f}%{'':8s} {s['coverage_near']:.0f}%")

    total_all = sum(s["total"] for s in summary.values())
    exact_all = sum(s["exact"] for s in summary.values())
    near_all = sum(s["exact"] + s["near"] for s in summary.values())
    print(f"\n{'TOTAL':<15} {total_all:<8} {exact_all:<8} {near_all:<8} {exact_all/total_all*100:.0f}%{'':8s} {near_all/total_all*100:.0f}%")

    return summary, results


if __name__ == "__main__":
    run_conservation_analysis()
