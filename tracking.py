#!/usr/bin/env python3
"""
Sistema de tracking de progresso do SmartLab BacEnd.
Rastreia quais genes ja passaram por cada etapa do pipeline.
"""

import os
import csv
from datetime import datetime

from config import TARGETS, SEQUENCES_DIR, GUIDES_DIR, PRIMERS_DIR, REPORTS_DIR, TRACKING_CSV


STEPS = ["fetch", "guides", "primers", "blast", "panel"]


def init_tracking():
    """Inicializa ou atualiza o CSV de tracking com todos os targets."""
    os.makedirs(REPORTS_DIR, exist_ok=True)
    existing = load_tracking()

    rows = []
    for gene_name in TARGETS:
        if gene_name in existing:
            rows.append(existing[gene_name])
        else:
            rows.append({
                "gene": gene_name,
                "fetch": "pending",
                "guides": "pending",
                "primers": "pending",
                "blast": "pending",
                "panel": "pending",
                "status": "queued",
                "updated": "",
            })

    save_tracking(rows)
    return rows


def load_tracking() -> dict:
    """Carrega tracking CSV como dict {gene: row}."""
    if not os.path.exists(TRACKING_CSV):
        return {}
    tracking = {}
    with open(TRACKING_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            tracking[row["gene"]] = row
    return tracking


def save_tracking(rows: list):
    """Salva tracking CSV."""
    fieldnames = ["gene", "fetch", "guides", "primers", "blast", "panel", "status", "updated"]
    with open(TRACKING_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def update_step(gene_name: str, step: str, result: str = "done"):
    """Atualiza o status de um step para um gene."""
    tracking = load_tracking()
    if gene_name not in tracking:
        tracking[gene_name] = {
            "gene": gene_name, "fetch": "pending", "guides": "pending",
            "primers": "pending", "blast": "pending", "panel": "pending",
            "status": "queued", "updated": "",
        }

    tracking[gene_name][step] = result
    tracking[gene_name]["updated"] = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Atualizar status geral
    steps_status = [tracking[gene_name].get(s, "pending") for s in STEPS]
    if all(s == "done" for s in steps_status):
        tracking[gene_name]["status"] = "complete"
    elif any(s == "failed" for s in steps_status):
        tracking[gene_name]["status"] = "failed"
    elif any(s == "done" for s in steps_status):
        tracking[gene_name]["status"] = "in_progress"
    else:
        tracking[gene_name]["status"] = "queued"

    save_tracking(list(tracking.values()))


def auto_detect_status():
    """Detecta automaticamente o status de cada gene baseado nos arquivos existentes."""
    tracking = load_tracking()
    if not tracking:
        tracking = {g: {"gene": g, "fetch": "pending", "guides": "pending",
                        "primers": "pending", "blast": "pending", "panel": "pending",
                        "status": "queued", "updated": ""} for g in TARGETS}

    for gene_name in TARGETS:
        if gene_name not in tracking:
            tracking[gene_name] = {
                "gene": gene_name, "fetch": "pending", "guides": "pending",
                "primers": "pending", "blast": "pending", "panel": "pending",
                "status": "queued", "updated": "",
            }

        row = tracking[gene_name]

        # Check fetch
        fasta_path = os.path.join(SEQUENCES_DIR, f"{gene_name}.fasta")
        row["fetch"] = "done" if os.path.exists(fasta_path) else "pending"

        # Check guides
        guides_path = os.path.join(GUIDES_DIR, f"{gene_name}_cas12a_guides.tsv")
        row["guides"] = "done" if os.path.exists(guides_path) else "pending"

        # Check primers
        primers_path = os.path.join(PRIMERS_DIR, f"{gene_name}_rpa_primers.tsv")
        row["primers"] = "done" if os.path.exists(primers_path) else "pending"

        # Check blast (look in specificity report)
        spec_path = os.path.join(REPORTS_DIR, "specificity_report.tsv")
        if os.path.exists(spec_path):
            with open(spec_path) as f:
                content = f.read()
                row["blast"] = "done" if gene_name in content else "pending"
        else:
            row["blast"] = "pending"

        # Check panel
        row["panel"] = "done" if row["fetch"] == "done" and row["guides"] == "done" and row["primers"] == "done" else "pending"

        # Status geral
        steps = [row[s] for s in STEPS]
        if all(s == "done" for s in steps):
            row["status"] = "complete"
        elif any(s == "done" for s in steps):
            row["status"] = "in_progress"
        else:
            row["status"] = "queued"

        row["updated"] = datetime.now().strftime("%Y-%m-%d %H:%M")

    save_tracking(list(tracking.values()))
    return tracking


def get_queued_genes(n: int = 5) -> list:
    """Retorna os proximos N genes com status 'queued'."""
    tracking = load_tracking()
    queued = [g for g, row in tracking.items() if row.get("status") == "queued"]
    return queued[:n]


def get_incomplete_genes() -> list:
    """Retorna genes que nao estao completos."""
    tracking = load_tracking()
    return [g for g, row in tracking.items() if row.get("status") != "complete"]


def print_status():
    """Imprime status atual de todos os genes."""
    tracking = load_tracking()
    if not tracking:
        print("Nenhum tracking encontrado. Execute init_tracking() primeiro.")
        return

    print(f"\n{'='*80}")
    print(f"SmartLab BacEnd - Status de Sequenciamento")
    print(f"{'='*80}")

    complete = sum(1 for r in tracking.values() if r.get("status") == "complete")
    in_prog = sum(1 for r in tracking.values() if r.get("status") == "in_progress")
    queued = sum(1 for r in tracking.values() if r.get("status") == "queued")
    total = len(tracking)

    print(f"\nResumo: {complete}/{total} completos | {in_prog} em progresso | {queued} na fila\n")

    print(f"{'Gene':<15} {'Fetch':<8} {'Guides':<8} {'Primers':<8} {'BLAST':<8} {'Panel':<8} {'Status':<12} {'Atualizado':<16}")
    print("-" * 95)

    for gene, row in sorted(tracking.items(), key=lambda x: x[1].get("status", "z")):
        status_icon = {"complete": "✓", "in_progress": "→", "queued": "○", "failed": "✗"}.get(row.get("status", ""), "?")
        print(
            f"{gene:<15} "
            f"{row.get('fetch','?'):<8} "
            f"{row.get('guides','?'):<8} "
            f"{row.get('primers','?'):<8} "
            f"{row.get('blast','?'):<8} "
            f"{row.get('panel','?'):<8} "
            f"{status_icon} {row.get('status',''):<10} "
            f"{row.get('updated',''):<16}"
        )


def main():
    print("Inicializando tracking...")
    init_tracking()
    print("Detectando status automaticamente...")
    auto_detect_status()
    print_status()


if __name__ == "__main__":
    main()
