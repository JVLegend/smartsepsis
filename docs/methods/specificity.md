# Specificity

!!! warning "Computational design phase"
    Specificity here means **in silico specificity**. Experimental
    cross-reactivity has not been measured.

## In silico approach

Candidate guides and primer pairs are screened against public reference
repositories using BLAST (`specificity_check.py`). Matches are filtered by
identity, coverage, and taxonomy to flag potential off-target hits.

## Planned stress-test

A planned stress-test pass will probe selectivity against:

- Close non-target genomes within the same genus.
- Common ocular surface commensals.
- Phylogenetically adjacent resistance gene variants.

Details will live at `scripts/README_stress_test.md` once the stress-test
suite lands.[^stress]

[^stress]: At the time of writing, `scripts/README_stress_test.md` is being
    prepared in a sibling work-stream and may not yet be present on `main`.

## Caveats

- BLAST-based specificity is a **necessary but not sufficient** check.
  Wet-lab cross-reactivity testing is required before any performance claim.
- The reference repositories evolve; specificity calls are valid as of the
  pull date recorded in pipeline outputs.

## Source

- `specificity_check.py` — BLAST orchestration and filtering.
