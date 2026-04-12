#!/usr/bin/env python3
"""
Orquestrador de batch do SmartLab BacEnd.
Roda o pipeline completo para N genes por vez.
Uso: python run_batch.py [N] [--priority P1,P2]
"""

import sys
import os
import time

# Adicionar o diretorio do projeto ao path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import load_targets, TARGETS, SEQUENCES_DIR, GUIDES_DIR, PRIMERS_DIR, REPORTS_DIR
from tracking import init_tracking, auto_detect_status, get_queued_genes, update_step, print_status


def run_fetch_for_gene(gene_name: str, target: dict) -> bool:
    """Roda fetch_sequences para um gene."""
    from fetch_sequences import fetch_gene_sequence, save_fasta
    accession = target.get("gene_accession", "")
    if not accession:
        print(f"  ✗ {gene_name}: sem gene_accession")
        return False

    fasta_path = os.path.join(SEQUENCES_DIR, f"{gene_name}.fasta")
    if os.path.exists(fasta_path):
        print(f"  ✓ {gene_name}: ja baixado")
        return True

    os.makedirs(SEQUENCES_DIR, exist_ok=True)
    content = fetch_gene_sequence(accession, gene_name)
    if content:
        save_fasta(content, fasta_path)
        return True
    return False


def run_guides_for_gene(gene_name: str) -> bool:
    """Roda design_guides para um gene."""
    from design_guides import design_guides_for_gene, save_guides
    from utils import parse_fasta

    fasta_path = os.path.join(SEQUENCES_DIR, f"{gene_name}.fasta")
    if not os.path.exists(fasta_path):
        return False

    os.makedirs(GUIDES_DIR, exist_ok=True)
    seqs = parse_fasta(fasta_path)
    if not seqs:
        return False

    sequence = list(seqs.values())[0]
    guides = design_guides_for_gene(gene_name, sequence)
    if guides:
        save_guides(gene_name, guides)
        return True
    return False


def run_primers_for_gene(gene_name: str) -> bool:
    """Roda design_primers para um gene."""
    from design_primers import load_best_guide, design_primers_for_target, save_primers
    from utils import parse_fasta

    guide = load_best_guide(gene_name)
    if not guide:
        return False

    fasta_path = os.path.join(SEQUENCES_DIR, f"{gene_name}.fasta")
    if not os.path.exists(fasta_path):
        return False

    os.makedirs(PRIMERS_DIR, exist_ok=True)
    seqs = parse_fasta(fasta_path)
    sequence = list(seqs.values())[0]
    guide_pos = int(guide["position"])

    primers = design_primers_for_target(gene_name, sequence, guide_pos)
    if primers:
        save_primers(gene_name, primers)
        return True
    return False


def run_pipeline_for_genes(gene_names: list, skip_blast: bool = False):
    """Roda o pipeline completo para uma lista de genes."""
    targets = TARGETS

    print(f"\n{'='*60}")
    print(f"SmartLab BacEnd - Batch Pipeline")
    print(f"Genes: {', '.join(gene_names)}")
    print(f"{'='*60}")

    for gene_name in gene_names:
        target = targets.get(gene_name, {})
        if not target:
            print(f"\n✗ {gene_name}: nao encontrado em TARGETS")
            continue

        print(f"\n{'─'*40}")
        print(f"[{gene_name}] ({target.get('pathogen', '')})")

        # Step 1: Fetch
        print(f"  1/4 Fetch...")
        ok = run_fetch_for_gene(gene_name, target)
        update_step(gene_name, "fetch", "done" if ok else "failed")
        if not ok:
            print(f"  ✗ Fetch falhou, pulando gene")
            continue
        time.sleep(0.5)

        # Step 2: Guides
        print(f"  2/4 Design guides...")
        ok = run_guides_for_gene(gene_name)
        update_step(gene_name, "guides", "done" if ok else "failed")
        if not ok:
            print(f"  ✗ Guides falhou")
            continue

        # Step 3: Primers
        print(f"  3/4 Design primers...")
        ok = run_primers_for_gene(gene_name)
        update_step(gene_name, "primers", "done" if ok else "failed")
        if not ok:
            print(f"  ✗ Primers falhou")
            continue

        # Step 4: BLAST (skip if requested)
        if skip_blast:
            print(f"  4/4 BLAST: pulado (--skip-blast)")
            update_step(gene_name, "blast", "pending")
        else:
            print(f"  4/4 BLAST: sera executado em batch ao final")
            update_step(gene_name, "blast", "pending")

        update_step(gene_name, "panel", "done")
        print(f"  ✓ {gene_name} processado!")

    print_status()


def main():
    # Parse args
    n = 4  # default batch size
    skip_blast = "--skip-blast" in sys.argv
    priority = None

    for arg in sys.argv[1:]:
        if arg.isdigit():
            n = int(arg)
        elif arg.startswith("--priority="):
            priority = arg.split("=")[1].split(",")

    # Init tracking
    init_tracking()
    auto_detect_status()

    # Get queued genes
    queued = get_queued_genes(n)
    if not queued:
        print("Nenhum gene na fila. Todos ja foram processados!")
        print_status()
        return

    print(f"Proximos {len(queued)} genes na fila: {queued}")
    run_pipeline_for_genes(queued, skip_blast=skip_blast)


if __name__ == "__main__":
    main()
