# Changelog

All notable changes to **SmartSepsis-Oph** are documented here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0-phase0] - 2026-05-11

### Added
- Initial public positioning as a **Phase 0 AI-first computational design
  program** for CRISPR-based ophthalmic point-of-care diagnostics.
- Repository metadata: `pyproject.toml` (PEP 621), `LICENSE` (MIT),
  `CONTRIBUTING.md`, `CITATION.cff`, `requirements.txt`,
  `requirements-dev.txt`, `MIGRATION.md`, `CHANGELOG.md`.
- Smoke test suite under `tests/` validating `pyproject.toml` parseability
  and import-health of root-level pipeline modules.
- Public website (`public/`) re-aligned with the computational-design phase
  scope (no clinical claims; design-stage language only).
- HuggingFace dataset published as **`smartsepsis-oph`** (panel + extended
  configurations) covering target sequences, designed Cas12a guides, RPA
  primers, and resistance-gene annotations.
- Pipeline modules (script-style, repo root):
  - **Sequence ingest & QC**: `fetch_sequences.py`, `conservation_analysis.py`
  - **Guide / probe design**: `design_guides.py`, `covariance_probes.py`,
    `crrna_secondary_structure.py`, `specificity_check.py`
  - **Primer design**: `design_primers.py`
  - **Multiplex assembly**: `multiplex_panel.py`
  - **Scoring (functional / language-model)**: `evo2_scoring.py`,
    `protein_scoring.py`, `protein_distance_matrix.py`,
    `prott5_ensemble.py`, `amrfinderplus_embed.py`,
    `phenotype_probe.py`, `phenotype_probe_v2.py`
  - **Structure prediction & features**: `protein_structure.py`,
    `structure_features_v3.py`, `nanofold_benchmark.py`,
    `nanofold_calibration_v5.py`
  - **Clinical interpretation / knowledge integration**:
    `card_integration.py`, `clinical_interpreter.py`,
    `tier2_pubmed_mining.py`
  - **Dataset build**: `build_hf_dataset.py`, `build_hf_dataset_extended.py`,
    `prodigal_to_prokka_gff.py`
  - **Orchestration & utilities**: `run_batch.py`, `tracking.py`,
    `config.py`, `utils.py`

### Notes
- The repository remains in **script-style layout** at the root for Phase 0.
  Migration to `src/smartsepsis_oph/…` is scheduled post-preprint — see
  [`MIGRATION.md`](./MIGRATION.md).
- Experimental (wet-lab) validation is **pending**; all current outputs are
  in silico designs.

[0.1.0-phase0]: https://github.com/JVLegend/smartsepsis/releases/tag/v0.1.0-phase0
