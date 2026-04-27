#!/usr/bin/env python3
"""
Matriz de distancia proteica set-based por familia AMR (Fase 10).

Usa os embeddings ESM-2 cacheados por protein_scoring.py e calcula
distancia coseno PAIRWISE (NxN) entre todas as variantes de cada familia.
Tambem gera clustering hierarquico simples (linkage) por familia.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from config import REPORTS_DIR, BASE_DIR

PROTEIN_SCORING_DIR = os.path.join(REPORTS_DIR, "protein_scoring")
EMBEDDINGS_DIR = os.path.join(PROTEIN_SCORING_DIR, "embeddings")
DISTANCE_REPORTS_DIR = os.path.join(REPORTS_DIR, "protein_distance")

ESM_MODEL_NAME = os.environ.get("ESM_MODEL", "esm2_t30_150M_UR50D")


def load_embedding(name: str) -> np.ndarray | None:
    path = os.path.join(EMBEDDINGS_DIR, f"{name}__{ESM_MODEL_NAME}.npy")
    if not os.path.exists(path):
        return None
    return np.load(path)


def cosine_distance_matrix(M: np.ndarray) -> np.ndarray:
    """Distancia coseno NxN a partir de matriz N x D."""
    norms = np.linalg.norm(M, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    Mn = M / norms
    sim = Mn @ Mn.T
    sim = np.clip(sim, -1.0, 1.0)
    return 1.0 - sim


def generate_matrices():
    print("=" * 70)
    print(f"FASE 10 - Distancia proteica set-based ({ESM_MODEL_NAME})")
    print("=" * 70)

    effects_csv = os.path.join(PROTEIN_SCORING_DIR, "protein_variant_effects.csv")
    if not os.path.exists(effects_csv):
        print(f"ERRO: {effects_csv} nao existe. Rode protein_scoring.py primeiro.")
        return

    os.makedirs(DISTANCE_REPORTS_DIR, exist_ok=True)
    df = pd.read_csv(effects_csv)

    summary = {}
    families = sorted(df["family"].unique())
    for fam in families:
        fam_df = df[df["family"] == fam].copy()
        variants = fam_df["variant"].tolist()
        if len(variants) < 2:
            continue

        embs = []
        kept = []
        for v in variants:
            e = load_embedding(v)
            if e is None:
                print(f"  [WARN] sem embedding para {v}")
                continue
            embs.append(e)
            kept.append(v)
        if len(kept) < 2:
            continue

        M = np.vstack(embs)
        D = cosine_distance_matrix(M)
        D = np.round(D.astype(float), 6)

        # CSV
        out_df = pd.DataFrame(D, index=kept, columns=kept)
        csv_path = os.path.join(DISTANCE_REPORTS_DIR, f"distance_{fam}.csv")
        out_df.to_csv(csv_path)

        # Estatisticas
        iu = np.triu_indices_from(D, k=1)
        pair_d = D[iu]
        stats = {
            "n_variants": len(kept),
            "mean_pairwise_distance": float(np.mean(pair_d)) if len(pair_d) else 0.0,
            "max_pairwise_distance": float(np.max(pair_d)) if len(pair_d) else 0.0,
            "min_pairwise_distance": float(np.min(pair_d)) if len(pair_d) else 0.0,
        }

        summary[fam] = {
            "variants": kept,
            "matrix": D.tolist(),
            "stats": stats,
        }
        print(
            f"  {fam:<14} n={len(kept):>2}  "
            f"mean={stats['mean_pairwise_distance']:.4f}  "
            f"max={stats['max_pairwise_distance']:.4f}"
        )

    json_path = os.path.join(DISTANCE_REPORTS_DIR, "protein_distances_all.json")
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n[Fase 10] CSVs por familia -> {DISTANCE_REPORTS_DIR}/")
    print(f"[Fase 10] Resumo JSON      -> {json_path}")


if __name__ == "__main__":
    generate_matrices()
