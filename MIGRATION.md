# Migration Plan — toward a `src/smartsepsis_oph/` package layout

## Current state (Phase 0)

The SmartSepsis-Oph repository currently keeps its Python modules at the **repo
root** as flat scripts (`design_guides.py`, `design_primers.py`, `run_batch.py`,
`evo2_scoring.py`, `protein_structure.py`, etc.). This reflects the project's
origin as a research prototype: each script is a self-contained pipeline stage
that imports a small `config.py` and a `utils.py` and can be invoked directly
via `python design_guides.py …`.

This layout is convenient for fast iteration, but it does **not** scale to a
distributable package, complicates testing, and prevents clean entry-point
registration in `pyproject.toml`.

## Why we are not refactoring yet

The Phase 0 deliverables prioritized for the next two quarters are:

1. **Preprint** describing the AI-first design program.
2. **Public dataset release** on HuggingFace (`smartsepsis-oph` — panel +
   extended).
3. **Partner / advisor outreach** (Mass Eye and Ear, HC-FMUSP collaborators,
   grant evaluators).

A package-layout refactor in the middle of this window would (a) churn import
paths across every notebook and downstream consumer, (b) invalidate
already-shared file references, and (c) cost reviewer-attention bandwidth with
no scientific benefit. The refactor is therefore **scheduled post-preprint**.

## Target layout

```
src/
  smartsepsis_oph/
    __init__.py
    config.py
    utils.py
    pipeline/
      __init__.py
      fetch_sequences.py
      design_guides.py
      covariance_probes.py
      design_primers.py
      specificity_check.py
      crrna_secondary_structure.py
      conservation_analysis.py
      run_batch.py
      tracking.py
    scoring/
      __init__.py
      evo2_scoring.py
      protein_scoring.py
      protein_distance_matrix.py
      prott5_ensemble.py
      amrfinderplus_embed.py
      phenotype_probe.py
      phenotype_probe_v2.py
    structure/
      __init__.py
      protein_structure.py
      structure_features_v3.py
      nanofold_benchmark.py
      nanofold_calibration_v5.py
    interpretation/
      __init__.py
      card_integration.py
      clinical_interpreter.py
      tier2_pubmed_mining.py
    multiplex/
      __init__.py
      multiplex_panel.py
    data/
      __init__.py
      build_hf_dataset.py
      build_hf_dataset_extended.py
      prodigal_to_prokka_gff.py
tests/
  ...
```

Each existing root `*.py` maps **one-to-one** into the tree above; no module
will be split or renamed during the move.

## Migration checklist (post-preprint)

- [ ] Create `src/smartsepsis_oph/` and the subpackage directories above, each
      with an `__init__.py`.
- [ ] `git mv` every root `*.py` into its target subpackage. Use `git mv` so
      Git preserves history.
- [ ] Rewrite intra-project imports:
      - `from config import …` → `from smartsepsis_oph.config import …`
      - `from utils import …` → `from smartsepsis_oph.utils import …`
      - `from tracking import …` → `from smartsepsis_oph.pipeline.tracking import …`
- [ ] Add `[project.scripts]` entry points in `pyproject.toml` for the
      argparse-driven mains (`design_guides`, `design_primers`, `run_batch`,
      `specificity_check`, etc.).
- [ ] Update `tests/test_smoke.py` to discover modules under
      `src/smartsepsis_oph/` instead of the repo root.
- [ ] Update CI workflow(s) to `pip install -e .[dev]` and run `pytest`.
- [ ] Update `README.md` usage examples to the new `python -m
      smartsepsis_oph.pipeline.design_guides …` invocation style.
- [ ] Update `docs/` (mkdocs) navigation and `mkdocstrings` import paths.
- [ ] Tag the release as `v0.2.0` once green.
