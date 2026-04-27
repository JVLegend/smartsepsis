#!/usr/bin/env python3
"""
Phenotype Probe (Fase 9) — supervised learning sobre embeddings ESM-2.

Treina classificadores multi-label (LogisticRegression OvR) prevendo
a partir do embedding proteico:
  - drug_class (multi-label: carbapenem, cephalosporin, penam, etc)
  - mechanism  (single-label: inactivation, target replacement, ...)

Avalia via Leave-One-Out CV (dataset pequeno: ~42 variantes / 12 familias).
Inspirado no workflow Learn da OpenProtein.AI.
"""

from __future__ import annotations

import csv
import json
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import MultiLabelBinarizer

from config import REPORTS_DIR

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

PROTEIN_SCORING_DIR = os.path.join(REPORTS_DIR, "protein_scoring")
EMBEDDINGS_DIR = os.path.join(PROTEIN_SCORING_DIR, "embeddings")
PHENOTYPE_REPORTS_DIR = os.path.join(REPORTS_DIR, "phenotype_probe")

ESM_MODEL_NAME = os.environ.get("ESM_MODEL", "esm2_t30_150M_UR50D")

# Labels CARD por familia (drug_class multi-label, mechanism single)
FAMILY_LABELS = {
    "mecA":       (["penam", "cephalosporin", "methicillin"], "antibiotic target replacement"),
    "mecA1":      (["penam", "cephalosporin", "methicillin"], "antibiotic target replacement"),
    "mecA2":      (["penam", "cephalosporin", "methicillin"], "antibiotic target replacement"),
    "blaKPC":     (["carbapenem", "cephalosporin", "penam"], "antibiotic inactivation"),
    "blaNDM":     (["carbapenem", "cephalosporin", "penam"], "antibiotic inactivation"),
    "blaOXA":     (["carbapenem", "cephalosporin", "penam"], "antibiotic inactivation"),
    "blaVIM":     (["carbapenem", "cephalosporin", "penam"], "antibiotic inactivation"),
    "blaIMP":     (["carbapenem", "cephalosporin", "penam"], "antibiotic inactivation"),
    "blaGES":     (["carbapenem", "cephalosporin"],          "antibiotic inactivation"),
    "blaCTX-M":   (["cephalosporin", "penam"],               "antibiotic inactivation"),
    "blaCTX-M-15":(["cephalosporin", "penam"],               "antibiotic inactivation"),
    "vanA":       (["glycopeptide"],                         "antibiotic target alteration"),
    "mcr-1":      (["polymyxin", "peptide"],                 "antibiotic target alteration"),
    "mcr-5":      (["polymyxin", "peptide"],                 "antibiotic target alteration"),
    "qnrS":       (["fluoroquinolone"],                      "antibiotic target protection"),
    "armA":       (["aminoglycoside"],                       "antibiotic target alteration"),
}


def load_embedding(name: str) -> np.ndarray | None:
    path = os.path.join(EMBEDDINGS_DIR, f"{name}__{ESM_MODEL_NAME}.npy")
    if not os.path.exists(path):
        return None
    return np.load(path)


def family_to_labels(fam: str) -> tuple[list[str], str] | None:
    if fam in FAMILY_LABELS:
        return FAMILY_LABELS[fam]
    # tenta prefixo
    for k, v in FAMILY_LABELS.items():
        if fam.startswith(k):
            return v
    return None


def run_phenotype_probe():
    print("=" * 70)
    print(f"FASE 9 - Phenotype probe (ESM-2 + sklearn) [{ESM_MODEL_NAME}]")
    print("=" * 70)

    os.makedirs(PHENOTYPE_REPORTS_DIR, exist_ok=True)

    effects_csv = os.path.join(PROTEIN_SCORING_DIR, "protein_variant_effects.csv")
    if not os.path.exists(effects_csv):
        print(f"ERRO: {effects_csv} nao existe. Rode protein_scoring.py primeiro.")
        return

    df = pd.read_csv(effects_csv)

    # Carregar embeddings + labels
    X, names, fams, drug_labels, mech_labels = [], [], [], [], []
    for _, row in df.iterrows():
        name = row["variant"]
        fam = row["family"]
        emb = load_embedding(name)
        if emb is None:
            continue
        labs = family_to_labels(fam)
        if labs is None:
            print(f"  [SKIP] sem labels CARD para familia {fam}")
            continue
        drugs, mech = labs
        X.append(emb)
        names.append(name)
        fams.append(fam)
        drug_labels.append(drugs)
        mech_labels.append(mech)

    if len(X) < 5:
        print(f"ERRO: amostras insuficientes ({len(X)}).")
        return

    X = np.vstack(X)
    print(f"Dataset: N={len(X)}  dim={X.shape[1]}  familias={len(set(fams))}")

    # === Drug class: multi-label OvR ===
    mlb = MultiLabelBinarizer()
    Y = mlb.fit_transform(drug_labels)
    print(f"Drug classes: {list(mlb.classes_)}")

    clf_drug = OneVsRestClassifier(LogisticRegression(max_iter=2000, C=1.0))
    clf_drug.fit(X, Y)
    Y_pred_train = clf_drug.predict(X)

    # LOO CV
    print("\nAvaliacao Leave-One-Out (drug class)...")
    loo = LeaveOneOut()
    Y_loo = np.zeros_like(Y)
    Y_loo_proba = np.zeros_like(Y, dtype=float)
    for train_idx, test_idx in loo.split(X):
        c = OneVsRestClassifier(LogisticRegression(max_iter=2000, C=1.0))
        c.fit(X[train_idx], Y[train_idx])
        Y_loo[test_idx] = c.predict(X[test_idx])
        try:
            Y_loo_proba[test_idx] = c.predict_proba(X[test_idx])
        except Exception:
            Y_loo_proba[test_idx] = Y_loo[test_idx]

    # Metricas
    exact_match = float(np.mean((Y_loo == Y).all(axis=1)))
    hamming = float(np.mean(Y_loo == Y))
    per_class_f1 = {}
    for j, cls in enumerate(mlb.classes_):
        tp = int(((Y_loo[:, j] == 1) & (Y[:, j] == 1)).sum())
        fp = int(((Y_loo[:, j] == 1) & (Y[:, j] == 0)).sum())
        fn = int(((Y_loo[:, j] == 0) & (Y[:, j] == 1)).sum())
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
        per_class_f1[cls] = round(f1, 3)
    print(f"  exact-match LOO = {exact_match:.3f}")
    print(f"  hamming acc LOO = {hamming:.3f}")
    print(f"  F1 per class    = {per_class_f1}")

    # === Mechanism: single-label ===
    mech_set = sorted(set(mech_labels))
    mech_idx = {m: i for i, m in enumerate(mech_set)}
    y_mech = np.array([mech_idx[m] for m in mech_labels])
    print(f"\nMechanisms: {mech_set}")

    # LOO so se houver >1 classe
    mech_acc_loo = None
    if len(mech_set) > 1:
        correct = 0
        for train_idx, test_idx in loo.split(X):
            if len(set(y_mech[train_idx])) < 2:
                continue
            c = LogisticRegression(max_iter=2000, C=1.0)
            c.fit(X[train_idx], y_mech[train_idx])
            pred = c.predict(X[test_idx])
            if pred[0] == y_mech[test_idx][0]:
                correct += 1
        mech_acc_loo = correct / len(X)
        print(f"  accuracy LOO = {mech_acc_loo:.3f}")

    clf_mech = LogisticRegression(max_iter=2000, C=1.0)
    if len(mech_set) > 1:
        clf_mech.fit(X, y_mech)

    # === Predicoes finais (treinado em tudo) ===
    rows = []
    for i, name in enumerate(names):
        pred_drug_idx = np.where(Y_pred_train[i] == 1)[0]
        pred_drugs = [mlb.classes_[k] for k in pred_drug_idx]
        try:
            proba_drug = clf_drug.predict_proba(X[i : i + 1])[0]
        except Exception:
            proba_drug = Y_pred_train[i].astype(float)

        if len(mech_set) > 1:
            pred_mech_i = clf_mech.predict(X[i : i + 1])[0]
            pred_mech = mech_set[pred_mech_i]
            mech_conf = float(np.max(clf_mech.predict_proba(X[i : i + 1])[0]))
        else:
            pred_mech = mech_set[0]
            mech_conf = 1.0

        rows.append(
            {
                "variant": name,
                "family": fams[i],
                "true_drug_class": ";".join(drug_labels[i]),
                "predicted_drug_class": ";".join(pred_drugs) if pred_drugs else "(none)",
                "drug_confidence_max": round(float(np.max(proba_drug)), 3),
                "true_mechanism": mech_labels[i],
                "predicted_mechanism": pred_mech,
                "mechanism_confidence": round(mech_conf, 3),
                "model": ESM_MODEL_NAME,
            }
        )

    csv_path = os.path.join(PHENOTYPE_REPORTS_DIR, "phenotype_predictions.csv")
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    metrics = {
        "model": ESM_MODEL_NAME,
        "n_samples": int(len(X)),
        "drug_classes": list(mlb.classes_),
        "drug_exact_match_loo": round(exact_match, 4),
        "drug_hamming_acc_loo": round(hamming, 4),
        "drug_per_class_f1_loo": per_class_f1,
        "mechanisms": mech_set,
        "mechanism_acc_loo": round(mech_acc_loo, 4) if mech_acc_loo is not None else None,
    }
    metrics_path = os.path.join(PHENOTYPE_REPORTS_DIR, "phenotype_metrics.json")
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    json_path = os.path.join(PHENOTYPE_REPORTS_DIR, "phenotype_predictions.json")
    with open(json_path, "w") as f:
        json.dump(rows, f, indent=2)

    print(f"\n[Fase 9] CSV     -> {csv_path}")
    print(f"[Fase 9] JSON    -> {json_path}")
    print(f"[Fase 9] Metrics -> {metrics_path}")


if __name__ == "__main__":
    run_phenotype_probe()
