#!/usr/bin/env python3
"""
Item #1 — Benchmark externo com NanoFold-public.

Demonstra que nosso pipeline ESM-2 (AMR-trained) generaliza para proteinas
nao-AMR usando o split validation do NanoFold (1000 cadeias diversas).

Pipeline:
  1. Baixa nanofold-public via HF datasets (1.09 GB total, val=~100MB)
  2. Embed val sequences (40-256aa) com ESM-2 (~10min M2)
  3. Aplica nosso classifier multi-label de drug class (treinado em AMRFinderPlus)
     sobre essas proteinas nao-AMR. Espera-se predicao baixa (none / unknown).
  4. Calcula taxa de "false positive" como metrica de calibracao
  5. Salva em reports/nanofold_benchmark/

Uso:
    /Users/iaparamedicos/envs/dev/bin/python nanofold_benchmark.py
    NANO_LIMIT=100 ... (sub-amostra rapida)
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import torch
import esm
from datasets import load_dataset
import joblib

from config import REPORTS_DIR

OUT_DIR = os.path.join(REPORTS_DIR, "nanofold_benchmark")
ESM_MODEL_NAME = os.environ.get("ESM_MODEL", "esm2_t30_150M_UR50D")
LIMIT = int(os.environ.get("NANO_LIMIT", "1000"))
ESM_MAX_LEN = 1022

MODEL_LAYERS = {"esm2_t30_150M_UR50D": 30, "esm2_t33_650M_UR50D": 33}

# Mapa AA-id -> letra (NanoFold usa AlphaFold ordering)
AF_AA = "ARNDCQEGHILKMFPSTWYV"


def aa_ids_to_seq(ids: list[int]) -> str:
    return "".join(AF_AA[i] if i < 20 else "X" for i in ids)


def load_esm():
    print(f"[ESM-2] {ESM_MODEL_NAME}")
    model, alphabet = esm.pretrained.load_model_and_alphabet(ESM_MODEL_NAME)
    model.eval()
    device = torch.device("cpu")
    model.to(device)
    return model, alphabet, MODEL_LAYERS[ESM_MODEL_NAME], device


def embed_seq(seq, model, alphabet, layer, device):
    bc = alphabet.get_batch_converter()
    seq = seq[:ESM_MAX_LEN]
    _, _, tokens = bc([("p", seq)])
    tokens = tokens.to(device)
    with torch.no_grad():
        out = model(tokens, repr_layers=[layer], return_contacts=False)
    rep = out["representations"][layer][0, 1 : 1 + len(seq)]
    return rep.mean(0).cpu().numpy().astype(np.float32)


def main():
    print("=" * 70)
    print("ITEM #1 — Benchmark externo NanoFold-public")
    print("=" * 70)
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Carregar NanoFold validation
    print("\nBaixando nanofold-public (1.09 GB, primeira vez)...")
    ds = load_dataset("ChrisHayduk/nanofold-public", split="validation")
    print(f"NanoFold val: {len(ds)} cadeias")

    # 2. Sub-amostrar
    if LIMIT > 0 and LIMIT < len(ds):
        ds = ds.shuffle(seed=42).select(range(LIMIT))
        print(f"Sub-amostra LIMIT={LIMIT}: {len(ds)} cadeias")

    # 3. Embed
    model, alphabet, layer, device = load_esm()
    print(f"\nEmbedding {len(ds)} sequencias com ESM-2 (~{len(ds) * 0.5:.0f}s estimado)...")

    X, names, lengths = [], [], []
    t0 = time.time()
    for i, row in enumerate(ds):
        seq = aa_ids_to_seq(row["aatype"])
        if len(seq) < 30:
            continue
        emb = embed_seq(seq, model, alphabet, layer, device)
        X.append(emb)
        names.append(row["chain_id"])
        lengths.append(len(seq))
        if (i + 1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  [{i + 1:>4}/{len(ds)}] elapsed={elapsed:.0f}s")
    print(f"[ESM-2] {len(X)} embeddings em {time.time() - t0:.0f}s")
    X = np.vstack(X)

    # 4. Aplicar classificador AMRFinderPlus (Fase 5a) — espera-se predicao "none"
    clf_path = os.path.join(REPORTS_DIR, "amrfinderplus", f"classifier_{ESM_MODEL_NAME}.joblib")
    if not os.path.exists(clf_path):
        print(f"AVISO: classifier nao encontrado em {clf_path}")
        print("       rode amrfinderplus_embed.py antes")
        return

    bundle = joblib.load(clf_path)
    clf, mlb = bundle["clf"], bundle["mlb"]
    print(f"\nClassifier: {len(mlb.classes_)} drug classes")

    pred = clf.predict(X)
    proba = clf.predict_proba(X)

    n_called = (pred.sum(axis=1) > 0).sum()
    avg_max_proba = float(proba.max(axis=1).mean())
    n_high_conf = int((proba.max(axis=1) > 0.5).sum())

    # Distribuicao de classes preditas (espera-se baixa)
    class_counts = pred.sum(axis=0)
    top5_idx = np.argsort(-class_counts)[:5]
    top5_classes = [(mlb.classes_[i], int(class_counts[i])) for i in top5_idx]

    metrics = {
        "n_chains": int(len(X)),
        "model": ESM_MODEL_NAME,
        "classifier": "AMRFinderPlus OvR LR",
        "false_positive_rate": round(float(n_called / len(X)), 4),
        "high_confidence_calls": n_high_conf,
        "avg_max_proba": round(avg_max_proba, 4),
        "top5_predicted_classes": top5_classes,
        "interpretation": (
            "Baixo false-positive rate = classifier nao confunde proteinas "
            "nao-AMR com AMR. Alto avg_proba indica falta de calibracao."
        ),
    }
    print(f"\n=== Resultados ===")
    print(f"  Cadeias avaliadas:     {metrics['n_chains']}")
    print(f"  False positive rate:   {metrics['false_positive_rate']:.2%}")
    print(f"  High confidence calls: {metrics['high_confidence_calls']}")
    print(f"  Avg max proba:         {metrics['avg_max_proba']:.3f}")
    print(f"  Top 5 false positives: {metrics['top5_predicted_classes']}")

    out_json = os.path.join(OUT_DIR, "metrics.json")
    with open(out_json, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\n[Saved] {out_json}")

    # Salvar embeddings p/ uso em item #2 (negativos)
    emb_path = os.path.join(OUT_DIR, f"embeddings_{ESM_MODEL_NAME}.npz")
    np.savez_compressed(
        emb_path,
        X=X,
        names=np.array(names, dtype=object),
        lengths=np.array(lengths, dtype=np.int32),
    )
    print(f"[Saved] {emb_path}  shape={X.shape}")
    print("\nDone. Use estes embeddings como negativos no item #2.")


if __name__ == "__main__":
    main()
