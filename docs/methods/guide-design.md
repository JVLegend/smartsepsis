# Guide-RNA design

!!! warning "Computational design phase"
    Guide-RNA designs in this repository are **in silico candidates**. No
    Cas12a collateral-cleavage validation has been performed. See
    [Roadmap — Phase 1](../project/roadmap.md#phase-1--wet-lab-feasibility).

## Approach

Cas12a guide-RNAs are designed by scanning reference target sequences for the
canonical **TTTV PAM** (V = A/C/G) and extracting the downstream spacer.
Candidate guides are filtered on:

- PAM context.
- Spacer composition (GC range, homopolymer avoidance).
- Predicted on-target activity.
- Position within the target gene.

See the source script: `design_guides.py`.

## Re-scoring with biophysical features

Candidate guides are then re-scored using **18 biophysical / sequence-context
features** via `covariance_probes.py`. The features capture local sequence
composition, predicted structure, and proxies for off-target behavior, and
are intended as a richer ranking signal than raw on-target score alone.

Details of the feature set are documented under [scoring](scoring.md).

## Target families

Initial library covers 12 resistance-gene families relevant to ocular
infection (see the [published library](../data/library.md) for the full
table). First declared validation target: ***S. aureus* / mecA**.

## Limitations

- All scoring is in silico; predicted activity has not been benchmarked
  against Cas12a cleavage data on these specific targets.
- Off-target profiling is BLAST-based — see [specificity](specificity.md) for
  caveats and planned stress-tests.

## Source

- `design_guides.py` — guide scanning and on-target scoring.
- `covariance_probes.py` — 18-feature re-scoring.
