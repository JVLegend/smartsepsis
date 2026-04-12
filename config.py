"""
Configuracao central do pipeline SmartLab BacEnd.
CRISPR-Cas12a paper-based diagnostic - IA para Medicos.
Suporta N alvos via targets_brazil.csv.
"""

import os
import csv

# === Diretorios ===
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SEQUENCES_DIR = os.path.join(BASE_DIR, "sequences")
GUIDES_DIR = os.path.join(BASE_DIR, "guides")
PRIMERS_DIR = os.path.join(BASE_DIR, "primers")
REPORTS_DIR = os.path.join(BASE_DIR, "reports")
TARGETS_CSV = os.path.join(BASE_DIR, "targets_brazil.csv")
TRACKING_CSV = os.path.join(BASE_DIR, "reports", "tracking_status.csv")


# === Alvos do painel (carregados do CSV) ===
def load_targets(csv_path: str = None, priority_filter: list = None) -> dict:
    """Carrega alvos do CSV. Opcionalmente filtra por prioridade (P1, P2, P3)."""
    path = csv_path or TARGETS_CSV
    if not os.path.exists(path):
        print(f"AVISO: {path} nao encontrado. Usando alvos padrao (mecA + blaKPC).")
        return {
            "mecA": {
                "name": "mecA",
                "organism": "Staphylococcus aureus",
                "gene_accession": "NG_047945.1",
                "pathogen": "MRSA",
                "clinical_relevance": "Methicillin-resistant S. aureus - major HAI pathogen",
                "priority": "P1",
            },
            "blaKPC": {
                "name": "blaKPC",
                "organism": "Klebsiella pneumoniae",
                "gene_accession": "NG_049243.1",
                "pathogen": "CRE",
                "clinical_relevance": "Carbapenem-resistant - WHO critical priority",
                "priority": "P1",
            },
        }

    targets = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if priority_filter and row.get("priority", "") not in priority_filter:
                continue
            targets[row["name"]] = {k: v for k, v in row.items()}
    return targets


TARGETS = load_targets()

# === Parâmetros Cas12a (LbCas12a / AsCas12a) ===
CAS12A = {
    "pam": "TTTV",  # V = A, C, G (não T)
    "pam_strand": "non-target",  # PAM reconhecido na fita não-alvo, upstream
    "spacer_length_min": 20,
    "spacer_length_max": 24,
    "spacer_length_optimal": 20,
    "gc_min": 0.30,
    "gc_max": 0.70,
    "gc_optimal_min": 0.40,
    "gc_optimal_max": 0.60,
    "max_homopolymer": 4,  # máximo de bases idênticas consecutivas
    "top_guides": 5,  # número de guides a reportar por gene
}

# === Parâmetros RPA (Recombinase Polymerase Amplification) ===
RPA = {
    "primer_length_min": 30,
    "primer_length_max": 35,
    "primer_length_optimal": 32,
    "amplicon_min": 100,
    "amplicon_max": 200,
    "amplicon_optimal": 150,
    "tm_min": 54.0,
    "tm_max": 67.0,
    "gc_min": 0.30,
    "gc_max": 0.70,
    "temperature": 37,  # °C
    "time_minutes": 20,
}

# === Reporter Cas12a ===
REPORTER = {
    "type": "ssDNA-FQ",
    "fluorophore": "FAM",
    "quencher": "BHQ-1",
    "sequence": "/56-FAM/TTATT/3BHQ_1/",  # 5nt ssDNA reporter typical
    "description": "ssDNA reporter clivado por atividade trans-cleavage de Cas12a",
}

# === Controles ===
CONTROLS = {
    "positive": {
        "name": "16S rRNA",
        "description": "Universal bacterial 16S rRNA - confirma presença de DNA bacteriano",
        "spot": "P",
    },
    "negative": {
        "name": "NTC",
        "description": "No-template control - confirma ausência de contaminação",
        "spot": "N",
    },
}

# === NCBI Entrez ===
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = "hackathon@example.com"  # required by NCBI policy

# === Cas12a Direct Repeat (scaffold para crRNA) ===
CAS12A_DIRECT_REPEAT = {
    "LbCas12a": "AATTTCTACTAAGTGTAGAT",
    "AsCas12a": "AATTTCTACTCTTGTAGAT",
}
