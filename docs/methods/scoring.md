# Scoring

!!! warning "Computational design phase"
    Scores are **predictive features**, not measured activities.

## Functional variant scoring

Variant-level functional impact is scored with an approach **inspired by
Evo 2 / EVEE**: large-scale sequence likelihood under a genomic language model
is used as a proxy for the fitness effect of a variant relative to wild type.

See `evo2_scoring.py`.

## 18 biophysical features (covariance probes)

Candidate guides are re-scored against 18 biophysical and sequence-context
features capturing:

- Local sequence composition (GC, dinucleotide bias).
- Predicted secondary structure context around the spacer.
- Position-weighted nucleotide identity.
- Proxies for off-target propensity.

These features are computed in `covariance_probes.py` and used to rank
candidates from [guide design](guide-design.md).

## CARD enrichment

Scored variants are cross-referenced against the **CARD** ontology
(`card_integration.py`) for:

- Drug class.
- Resistance mechanism.
- Ontology lineage.

This feeds the [clinical interpretation](pipeline-overview.md#7-clinical-interpretation--multiplex)
stage and the [published library](../data/library.md).

## Caveats

- Sequence-likelihood scores are **proxies**, not measured fitness.
- The 18-feature panel was built to rank within a candidate pool, not to give
  calibrated absolute activity.

## Sources

- `evo2_scoring.py` — functional variant scoring.
- `covariance_probes.py` — 18-feature re-scoring.
- `card_integration.py` — CARD ontology enrichment.
