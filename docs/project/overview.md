# Project overview

## What this project is

SmartSepsis-Oph is an AI-first **computational design program** for CRISPR-based
point-of-care diagnostics targeting ophthalmic infections — endophthalmitis,
microbial keratitis, and ocular surface infections.

The pipeline produces:

- CRISPR-Cas12a guide-RNA libraries for 12 resistance-gene families.
- Isothermal RPA primer designs compatible with thermocycler-free workflows.
- In silico specificity analyses against public reference repositories.
- Functional variant impact scoring (inspired by Evo 2 / EVEE).
- A paper-strip assay architecture under design.

First declared target: ***S. aureus* / mecA**.

## What this project is *not*

- **Not a medical device.**
- **Not validated experimentally** — every metric in this repository is in silico.
- **Not a substitute** for culture, antibiogram, or clonal surveillance.
- **Not approved** for any clinical use under ANVISA, FDA, or any other regulatory body.

Performance metrics (sensitivity, specificity, limit of detection, time-to-result)
will only be reported after experimental validation with a wet-lab partner.

## The ophthalmology gap

Existing molecular diagnostics are built for sample volumes orders of magnitude
larger than what ophthalmology yields:

- **Vitreous tap:** sub-microliter.
- **Corneal scrape / ocular swab:** trace amounts of nucleic acid.

Cultures take 48–72 hours; endophthalmitis blinds in 24–48. The diagnostic
window is incompatible with the standard of care.

## Why AI-first

We use protein language models (ESM-2, ProtT5), structural prediction
(ESMFold / ColabFold / AlphaFold), and functional variant scoring to design
guide-RNA and primer libraries that are evaluated computationally **before**
any wet-lab cycle is spent. This is the only realistic way to compress the
design-validation loop for a sample-scarce, time-critical indication.

## See also

- [Roadmap](roadmap.md) — phased plan.
- [Team](team.md) — who is involved.
- [Partners](partners.md) — what we are seeking.
- [Methods](../methods/pipeline-overview.md) — the seven-step pipeline.
- [Disclaimers](../disclaimers.md) — regulatory and research-use framing.
