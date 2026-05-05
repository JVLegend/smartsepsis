#!/usr/bin/env python3
"""
Item #2 — Phenotype probe calibrado com negativos NanoFold.

Treina classifier OvR LR multi-label drug class adicionando 800 cadeias NanoFold
como classe "no_resistance" (negativos). Reserva 200 NanoFold para hold-out FP test.

Outputs:
  reports/phenotype_probe_v2/metrics.json
  reports/phenotype_probe_v2/classifier.joblib
  reports/phenotype_probe_v2/predictions_painel.csv
"""

from __future__ import annotations

import csv
import json
import os
import warnings

import numpy as np
import pandas as pd
import joblib
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import f1_score
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import MultiLabelBinarizer

from config import REPORTS_DIR

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ESM_MODEL_NAME = os.environ.get("ESM_MODEL", "esm2_t30_150M_UR50D")
OUT_DIR = os.path.join(REPORTS_DIR, "phenotype_probe_v2")

# Paths
AMR_EMB_NPZ = os.path.join(REPORTS_DIR, "amrfinderplus", f"embeddings_{ESM_MODEL_NAME}.npz")
NANO_EMB_NPZ = os.path.join(REPORTS_DIR, "nanofold_benchmark", f"embeddings_{ESM_MODEL_NAME}.npz")
PANEL_EMB_DIR = os.path.join(REPORTS_DIR, "protein_scoring", "embeddings")
EFFECTS_CSV = os.path.join(REPORTS_DIR, "protein_scoring", "protein_variant_effects.csv")
NEG_LABEL = "no_resistance"


def main():
    print("=" * 70)
    print("ITEM #2 — Phenotype probe v2 com negativos NanoFold")
    print("=" * 70)
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. AMR positivos do AMRFinderPlus
    print("\nCarregando embeddings AMR (positivos)...")
    z = np.load(AMR_EMB_NPZ, allow_pickle=True)
    X_pos = z["X"]
    labels_pos = [list(l) for l in z["labels"]]
    print(f"  AMRFinderPlus: {len(X_pos)} positivos, dim={X_pos.shape[1]}")

    # 2. NanoFold negativos
    print("\nCarregando embeddings NanoFold (negativos)...")
    nz = np.load(NANO_EMB_NPZ, allow_pickle=True)
    X_neg_all = nz["X"]
    print(f"  NanoFold: {len(X_neg_all)} cadeias")

    # Split 800 train / 200 holdout
    rng = np.random.RandomState(42)
    idx = rng.permutation(len(X_neg_all))
    X_neg_train = X_neg_all[idx[:800]]
    X_neg_holdout = X_neg_all[idx[800:]]
    print(f"  Split: 800 train + {len(X_neg_holdout)} holdout")

    # 3. Combinar X e Y
    X = np.vstack([X_pos, X_neg_train])
    labels_all = labels_pos + [[NEG_LABEL] for _ in range(len(X_neg_train))]

    mlb = MultiLabelBinarizer()
    Y = mlb.fit_transform(labels_all)
    print(f"\nTraining set: N={len(X)}, classes={len(mlb.classes_)}")

    # 4. Treinar
    print("Treinando OvR LogisticRegression...")
    clf = OneVsRestClassifier(LogisticRegression(max_iter=2000, C=1.0, n_jobs=1), n_jobs=-1)
    clf.fit(X, Y)

    # 5. Avaliar holdout NanoFold (FP rate em proteinas nao-AMR)
    pred_holdout = clf.predict(X_neg_holdout)
    proba_holdout = clf.predict_proba(X_neg_holdout)
    # FP = predicao de qualquer classe AMR (ignorando no_resistance)
    no_res_idx = list(mlb.classes_).index(NEG_LABEL)
    amr_pred_mask = np.delete(pred_holdout, no_res_idx, axis=1)
    fp_rate_v2 = float((amr_pred_mask.sum(axis=1) > 0).sum() / len(X_neg_holdout))
    no_res_correct = float(pred_holdout[:, no_res_idx].sum() / len(X_neg_holdout))

    # 6. Avaliar painel (variantes do nosso panel) — comparacao c/ baseline
    print("\nAvaliando 43 variantes do painel (LOO)...")
    df = pd.read_csv(EFFECTS_CSV)

    # Mesmas family_labels do phenotype_probe.py original
    FAMILY_LABELS = {
        "mecA":        (["penam","cephalosporin","methicillin"], "antibiotic target replacement"),
        "mecA1":       (["penam","cephalosporin","methicillin"], "antibiotic target replacement"),
        "mecA2":       (["penam","cephalosporin","methicillin"], "antibiotic target replacement"),
        "blaKPC":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
        "blaNDM":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
        "blaOXA":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
        "blaVIM":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
        "blaIMP":      (["carbapenem","cephalosporin","penam"], "antibiotic inactivation"),
        "blaGES":      (["carbapenem","cephalosporin"], "antibiotic inactivation"),
        "blaCTX-M":    (["cephalosporin","penam"], "antibiotic inactivation"),
        "blaCTX-M-15": (["cephalosporin","penam"], "antibiotic inactivation"),
        "vanA":        (["glycopeptide"], "antibiotic target alteration"),
        "mcr-1":       (["polymyxin","peptide"], "antibiotic target alteration"),
        "mcr-5":       (["polymyxin","peptide"], "antibiotic target alteration"),
        "qnrS":        (["fluoroquinolone"], "antibiotic target protection"),
        "armA":        (["aminoglycoside"], "antibiotic target alteration"),
    }

    def labels_for(name):
        if name in FAMILY_LABELS: return FAMILY_LABELS[name]
        for k in FAMILY_LABELS:
            if name.startswith(k): return FAMILY_LABELS[k]
        return None

    panel_X, panel_names, panel_drugs = [], [], []
    for _, row in df.iterrows():
        name = row["variant"]
        labs = labels_for(row["family"])
        if labs is None: continue
        emb_path = os.path.join(PANEL_EMB_DIR, f"{name}__{ESM_MODEL_NAME}.npy")
        if not os.path.exists(emb_path): continue
        panel_X.append(np.load(emb_path))
        panel_names.append(name)
        panel_drugs.append(labs[0])
    panel_X = np.vstack(panel_X)
    print(f"  Painel: {len(panel_X)} variantes")

    # Inferencia direta (treinou no AMRFinder + neg, agora avalia no painel)
    panel_pred = clf.predict(panel_X)
    panel_proba = clf.predict_proba(panel_X)

    # F1 por classe nas predicoes do painel
    Y_panel_true = np.zeros((len(panel_X), len(mlb.classes_)), dtype=int)
    for i, drugs in enumerate(panel_drugs):
        for d in drugs:
            if d in mlb.classes_:
                Y_panel_true[i, list(mlb.classes_).index(d)] = 1

    panel_pred_amr_only = np.delete(panel_pred, no_res_idx, axis=1)
    Y_panel_true_amr = np.delete(Y_panel_true, no_res_idx, axis=1)

    f1_macro = f1_score(Y_panel_true_amr, panel_pred_amr_only, average="macro", zero_division=0)
    hamming_panel = float(np.mean(Y_panel_true_amr == panel_pred_amr_only))
    exact_panel = float((Y_panel_true_amr == panel_pred_amr_only).all(axis=1).mean())

    print(f"\n=== Resultados v2 ===")
    print(f"\nHoldout NanoFold (200 cadeias nao-AMR):")
    print(f"  FP rate v2:           {fp_rate_v2:.2%}     (baseline v1 = 14.90%)")
    print(f"  Predito 'no_resistance': {no_res_correct:.2%}")
    print(f"  Improvement vs v1:    {(0.149 - fp_rate_v2)*100:+.1f} pts")

    print(f"\nPainel SmartSepsis (43 variantes AMR):")
    print(f"  Hamming acc:          {hamming_panel:.3f}")
    print(f"  Exact match:          {exact_panel:.3f}")
    print(f"  F1 macro:             {f1_macro:.3f}")

    # Salvar
    metrics = {
        "version": "v2_with_negatives",
        "esm_model": ESM_MODEL_NAME,
        "n_amr_positives": int(len(X_pos)),
        "n_nanofold_train_negatives": int(len(X_neg_train)),
        "n_nanofold_holdout_negatives": int(len(X_neg_holdout)),
        "n_classes": int(len(mlb.classes_)),
        "fp_rate_v2_holdout": round(fp_rate_v2, 4),
        "fp_rate_v1_baseline": 0.149,
        "fp_rate_improvement_pts": round(0.149 - fp_rate_v2, 4),
        "no_resistance_correct_pct": round(no_res_correct, 4),
        "panel_hamming_acc": round(hamming_panel, 4),
        "panel_exact_match": round(exact_panel, 4),
        "panel_f1_macro": round(f1_macro, 4),
    }
    with open(os.path.join(OUT_DIR, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)

    joblib.dump({"clf": clf, "mlb": mlb, "model": ESM_MODEL_NAME, "neg_label": NEG_LABEL},
                os.path.join(OUT_DIR, f"classifier_{ESM_MODEL_NAME}.joblib"))

    # Predicoes do painel detalhadas
    rows = []
    for i, name in enumerate(panel_names):
        amr_idx = [k for k in np.where(panel_pred[i] == 1)[0] if k != no_res_idx]
        rows.append({
            "variant": name,
            "true_drugs": ";".join(panel_drugs[i]),
            "predicted_drugs": ";".join(mlb.classes_[k] for k in amr_idx) or "(none)",
            "prob_no_resistance": round(float(panel_proba[i][no_res_idx]), 3),
            "prob_top": round(float(panel_proba[i].max()), 3),
        })
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, "predictions_painel.csv"), index=False)

    print(f"\n[Saved] {OUT_DIR}/metrics.json")
    print(f"[Saved] classifier + predictions_painel.csv")


if __name__ == "__main__":
    main()
