#!/usr/bin/env python3
"""
Score funcional proteico de variantes AMR via ESM-2 (Meta AI).
Fase 7 — MELHORIAS_v2.md.

Pipeline:
  1. Para cada variante (DNA): pegar maior ORF (Table 11 bacteriana) -> proteina
  2. Embedding ESM-2 (mean-pooled, 640d com t30_150M)
  3. Cache em .npy (uma vez por sequencia)
  4. Distancia coseno entre variante e referencia da familia
  5. Classificacao funcional (synonymous / conserved / moderate / non_functional)

Hardware-aware: defaults para M2 16GB (CPU, modelo 150M).
Set ESM_MODEL=esm2_t33_650M_UR50D para usar 650M (precisa 32GB ideal).
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import requests
import torch
from Bio.Seq import Seq

import esm

from config import BASE_DIR, NCBI_BASE_URL, NCBI_EMAIL, REPORTS_DIR, SEQUENCES_DIR
from utils import parse_fasta, reverse_complement

# === Configuracao ===
PROTEIN_REPORTS_DIR = os.path.join(REPORTS_DIR, "protein_scoring")
EMBEDDINGS_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "embeddings")
PROTEINS_FASTA_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "proteins")

ESM_MODEL_NAME = os.environ.get("ESM_MODEL", "esm2_t30_150M_UR50D")
# Mapa modelo -> ultima camada (representacao)
MODEL_LAYERS = {
    "esm2_t6_8M_UR50D": 6,
    "esm2_t12_35M_UR50D": 12,
    "esm2_t30_150M_UR50D": 30,
    "esm2_t33_650M_UR50D": 33,
    "esm2_t36_3B_UR50D": 36,
}
ESM_MAX_LEN = 1022  # ESM-2 aceita ate 1024 com BOS/EOS


# ---------------------------------------------------------------------
# Tradução DNA -> proteína
# ---------------------------------------------------------------------
def longest_orf(dna: str, min_len: int = 50) -> str:
    """Retorna a maior ORF (proteina) considerando 6 frames, codigo bacteriano (Table 11).

    A ORF deve comecar em M e nao conter '*'. Retorna a proteina sem o stop final.
    """
    dna = dna.upper().replace("\n", "").replace(" ", "")
    best = ""
    for strand_seq in (dna, reverse_complement(dna)):
        for frame in range(3):
            sub = strand_seq[frame:]
            trim = len(sub) - (len(sub) % 3)
            if trim < 3:
                continue
            try:
                prot = str(Seq(sub[:trim]).translate(table=11))
            except Exception:
                continue
            for orf in prot.split("*"):
                idx = orf.find("M")
                if idx < 0:
                    continue
                candidate = orf[idx:]
                if len(candidate) >= min_len and len(candidate) > len(best):
                    best = candidate
    return best


# ---------------------------------------------------------------------
# Sequências NCBI
# ---------------------------------------------------------------------
def fetch_sequence_ncbi(accession: str) -> str | None:
    url = f"{NCBI_BASE_URL}/efetch.fcgi"
    params = {
        "db": "nucleotide",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "email": NCBI_EMAIL,
    }
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        content = resp.text.strip()
        if content.startswith(">"):
            lines = content.split("\n")
            return "".join(l.strip() for l in lines[1:] if not l.startswith(">")).upper()
    except Exception as e:
        print(f"    [NCBI] Erro fetch {accession}: {e}")
    return None


def ensure_dna(name: str, accession: str) -> str | None:
    fasta_path = os.path.join(SEQUENCES_DIR, f"{name}.fasta")
    if os.path.exists(fasta_path):
        seqs = parse_fasta(fasta_path)
        if seqs:
            return list(seqs.values())[0]
    print(f"  Baixando NCBI {name} ({accession})...")
    seq = fetch_sequence_ncbi(accession)
    if seq:
        os.makedirs(SEQUENCES_DIR, exist_ok=True)
        with open(fasta_path, "w") as f:
            f.write(f">{name} {accession}\n{seq}\n")
        time.sleep(0.4)
    return seq


# ---------------------------------------------------------------------
# ESM-2
# ---------------------------------------------------------------------
def load_esm_model():
    print(f"[ESM-2] Carregando modelo {ESM_MODEL_NAME} (download na 1a vez)...")
    model, alphabet = esm.pretrained.load_model_and_alphabet(ESM_MODEL_NAME)
    model.eval()
    # MPS tem bugs conhecidos com fair-esm; CPU é seguro e rápido o suficiente
    device = torch.device("cpu")
    model.to(device)
    layer = MODEL_LAYERS[ESM_MODEL_NAME]
    print(f"[ESM-2] Modelo pronto. Layer={layer}, device={device}")
    return model, alphabet, layer, device


def embed_protein(
    name: str,
    protein_seq: str,
    model,
    alphabet,
    layer: int,
    device,
) -> np.ndarray:
    """Mean-pooled embedding (sem BOS/EOS). Cacheia em .npy."""
    cache_path = os.path.join(EMBEDDINGS_DIR, f"{name}__{ESM_MODEL_NAME}.npy")
    if os.path.exists(cache_path):
        return np.load(cache_path)

    seq = protein_seq[:ESM_MAX_LEN]
    batch_converter = alphabet.get_batch_converter()
    _, _, tokens = batch_converter([(name, seq)])
    tokens = tokens.to(device)
    with torch.no_grad():
        out = model(tokens, repr_layers=[layer], return_contacts=False)
    rep = out["representations"][layer][0, 1 : 1 + len(seq)]  # remove BOS, e EOS
    emb = rep.mean(0).cpu().numpy().astype(np.float32)
    os.makedirs(EMBEDDINGS_DIR, exist_ok=True)
    np.save(cache_path, emb)
    return emb


def cosine_distance(a: np.ndarray, b: np.ndarray) -> float:
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na == 0 or nb == 0:
        return 1.0
    cos = float(np.dot(a, b) / (na * nb))
    cos = max(-1.0, min(1.0, cos))
    return 1.0 - cos


# ---------------------------------------------------------------------
# Classificação de impacto
# ---------------------------------------------------------------------
def classify_impact(prot_dist: float, dna_identical: bool, prot_identical: bool) -> tuple[str, str]:
    if dna_identical:
        return "IDENTICAL", "Sequencia identica a referencia"
    if prot_identical:
        return "SYNONYMOUS", "Mutacao sinonima (proteina identica)"
    if prot_dist < 0.01:
        return "CONSERVED", "Conservada (impacto minimo)"
    if prot_dist < 0.05:
        return "MISSENSE_LIGHT", "Missense leve"
    if prot_dist < 0.15:
        return "MISSENSE_MODERATE", "Missense moderado"
    return "NON_FUNCTIONAL", "Provavel perda de funcao / truncamento"


# ---------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------
def run_protein_scoring():
    print("=" * 70)
    print(f"FASE 7 - Score funcional proteico via ESM-2 ({ESM_MODEL_NAME})")
    print("SmartLab BacEnd | IA para Medicos")
    print("=" * 70)

    os.makedirs(PROTEIN_REPORTS_DIR, exist_ok=True)
    os.makedirs(EMBEDDINGS_DIR, exist_ok=True)
    os.makedirs(PROTEINS_FASTA_DIR, exist_ok=True)
    os.makedirs(SEQUENCES_DIR, exist_ok=True)

    targets_csv = os.path.join(BASE_DIR, "targets_brazil.csv")
    variants_csv = os.path.join(BASE_DIR, "targets_brazil_variants.csv")
    if not (os.path.exists(targets_csv) and os.path.exists(variants_csv)):
        print("ERRO: targets_brazil.csv ou targets_brazil_variants.csv ausente.")
        sys.exit(1)

    # Famílias de referência: do CSV de targets, prioridade nao-var
    family_ref = {}
    with open(targets_csv) as f:
        for row in csv.DictReader(f):
            if "var" in row.get("priority", "").lower():
                continue
            family_ref[row["name"]] = row["gene_accession"]

    # Variantes
    with open(variants_csv) as f:
        variants = list(csv.DictReader(f))

    # Garantir todas as DNAs em disco
    needed = {}
    for v in variants:
        needed[v["name"]] = v["gene_accession"]
        fam = v.get("gene_family") or v["name"]
        if fam in family_ref:
            needed[fam] = family_ref[fam]
        else:
            # familia pode coincidir com a propria variante de ref (ex: vanA)
            needed.setdefault(fam, v["gene_accession"])

    print(f"Garantindo {len(needed)} sequencias DNA em disco...")
    dna_seqs: dict[str, str] = {}
    for name, acc in needed.items():
        seq = ensure_dna(name, acc)
        if seq:
            dna_seqs[name] = seq
        else:
            print(f"  [WARN] sem DNA para {name}")

    # Traduzir tudo para proteína (longest ORF)
    print("\nTraduzindo DNA -> proteina (longest ORF, table 11)...")
    proteins: dict[str, str] = {}
    for name, dna in dna_seqs.items():
        prot = longest_orf(dna)
        if not prot:
            print(f"  [WARN] sem ORF para {name}")
            continue
        proteins[name] = prot
        out = os.path.join(PROTEINS_FASTA_DIR, f"{name}.fasta")
        with open(out, "w") as f:
            f.write(f">{name}\n{prot}\n")
    print(f"  {len(proteins)}/{len(dna_seqs)} traduzidas com sucesso")

    # ESM-2
    model, alphabet, layer, device = load_esm_model()

    print("\nGerando embeddings ESM-2 (cacheados em .npy)...")
    embeddings: dict[str, np.ndarray] = {}
    t0 = time.time()
    for i, (name, seq) in enumerate(proteins.items(), 1):
        emb = embed_protein(name, seq, model, alphabet, layer, device)
        embeddings[name] = emb
        print(f"  [{i:>2}/{len(proteins)}] {name:<18} L={len(seq):>4}aa dim={emb.shape[0]}")
    print(f"[ESM-2] {len(embeddings)} embeddings em {time.time() - t0:.1f}s")

    # Scoring
    print("\nCalculando distancia coseno variante <-> referencia da familia...")
    results = []
    for v in variants:
        name = v["name"]
        fam = v.get("gene_family") or name

        if name not in embeddings:
            continue
        if fam not in embeddings:
            # Sem referencia da familia -> usar a propria variante como ref (dist=0)
            ref_name = name
        else:
            ref_name = fam

        emb_var = embeddings[name]
        emb_ref = embeddings[ref_name]
        prot_var = proteins[name]
        prot_ref = proteins[ref_name]
        dna_var = dna_seqs[name]
        dna_ref = dna_seqs.get(ref_name, dna_var)

        prot_dist = cosine_distance(emb_var, emb_ref)
        impact, label = classify_impact(
            prot_dist,
            dna_identical=(dna_var == dna_ref),
            prot_identical=(prot_var == prot_ref),
        )

        results.append(
            {
                "variant": name,
                "family": fam,
                "ref_used": ref_name,
                "dna_length": len(dna_var),
                "protein_length": len(prot_var),
                "protein_distance_cosine": round(prot_dist, 6),
                "impact_type": impact,
                "label": label,
                "protein_score": round(prot_dist * 1000, 2),
                "model": ESM_MODEL_NAME,
            }
        )
        print(f"  {name:<18} fam={fam:<12} dist={prot_dist:.4f}  {impact}")

    # Persistencia
    if not results:
        print("Sem resultados.")
        return

    csv_path = os.path.join(PROTEIN_REPORTS_DIR, "protein_variant_effects.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)

    json_path = os.path.join(PROTEIN_REPORTS_DIR, "protein_scores.json")
    with open(json_path, "w") as f:
        json.dump({r["variant"]: r for r in results}, f, indent=2)

    print(f"\n[Fase 7] CSV   -> {csv_path}")
    print(f"[Fase 7] JSON  -> {json_path}")
    print(f"[Fase 7] Embeddings cache -> {EMBEDDINGS_DIR}/")


if __name__ == "__main__":
    run_protein_scoring()
