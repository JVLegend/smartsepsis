# Structural analysis

!!! warning "Computational design phase"
    Structures are **predicted**, not experimentally determined. Predicted
    structures carry inherent uncertainty (see pLDDT).

## Embeddings

Each variant in the [published library](../data/library.md) is annotated
with two protein-language-model embeddings:

- **ESM-2** (`esm2_t30_150M_UR50D`) — 640-dim, mean-pooled across residues.
- **ProtT5** (`Rostlab/prot_t5_xl_uniref50`) — 1024-dim, mean-pooled.

Concatenated, this gives a 1664-dim multimodal protein representation usable
as a downstream feature for classifiers and similarity search.

## 3D structure prediction

Per-variant 3D structures are predicted using a tiered approach:

| Method | Used for |
|---|---|
| **ESMFold** (HuggingFace transformers) | Proteins ≤ 400 aa. |
| **ColabFold AF2** | Mid-size targets (e.g. mcr-1, mcr-5 ~540 aa). |
| **AlphaFold Server (AF3)** | Larger targets (mecA1, mecA2 ~665 aa). |

`rank_1` is selected per variant.

## Structure descriptors

From the predicted Cα coordinates we compute:

- Radius of gyration (`struct_rg`).
- Compactness ratio `Rg / L^0.6` (`struct_compactness`).
- Contact density (fraction of Cα–Cα pairs < 8 Å).
- Mean Cα–Cα distance.
- Aspect ratio.
- **Mean pLDDT** — prediction confidence (0–100). Use this for downstream
  weighting; local regions can be substantially lower than the mean.

## Caveats

- AlphaFold Server (AF3) terms apply to the mecA1 / mecA2 PDBs — research
  use only, no commercial redistribution. For commercial use, regenerate
  via ColabFold AF2 or AlphaFold 2 / OpenFold.

## See also

- [Published library](../data/library.md) — schema and access.
- [Scoring](scoring.md) — how structural features feed ranking.
