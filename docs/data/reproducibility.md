# Reproducibility

!!! warning "Preliminary"
    Reproducibility documentation is **preliminary**. The pipeline scripts
    currently live at the repository root and are being migrated to a
    packaged `src/` layout ahead of the framework preprint. Expect rough
    edges until that migration lands. See `MIGRATION.md` in the repository
    for current status (file may not yet be present on `main`).

## Minimal flow

```bash
python fetch_sequences.py
python design_guides.py
python covariance_probes.py
python design_primers.py
python specificity_check.py
python run_batch.py 4
```

## What is reproducible today

- The **published library** ([data](library.md)) — fully reproducible from
  RefSeq accessions + pinned model checkpoints + the construction notes on
  Hugging Face.
- **Guide and primer designs** — deterministic given a reference sequence
  and the scripts at their pinned commit.
- **In silico specificity** — deterministic given a snapshot of the BLAST
  reference repositories on the pull date recorded in pipeline outputs.

## What is not yet reproducible

- **Performance metrics on real samples** — there are none. See
  [Disclaimers](../disclaimers.md).
- **Stress-test specificity** — planned, see
  [specificity](../methods/specificity.md).

## Asking for help

If you want to reproduce a specific result, open an issue at
<https://github.com/JVLegend/smartsepsis/issues> or email
**iaparamedicos@gmail.com**. We will help.

## See also

- [Pipeline overview](../methods/pipeline-overview.md).
- [Contributing](../contributing.md).
