# Disclaimers

!!! danger "Research use only"
    SmartSepsis-Oph is a research program by **IA para Médicos**. The code
    and designs in this repository are for **research use only**.

## No clinical claims

- **No claim of clinical performance is made.**
- All results in this repository are **in silico**.
- Performance metrics (sensitivity, specificity, limit of detection,
  time-to-result) will only be reported **after** experimental validation
  with a wet-lab partner.

## Not a medical device

- SmartSepsis-Oph is **not** a medical device.
- No medical device approval has been sought or granted under **ANVISA**,
  **FDA**, or any other regulatory body.
- Do **not** use any output of this pipeline to make clinical decisions.

## Regulatory status

- **ANVISA (Brazil):** no submission. The eventual pathway is anticipated
  under RDC 830/2023 (IVD).
- **FDA (USA):** no submission.
- **CE-IVDR (EU):** no submission.

Regulatory submission is a **Phase 4** activity on the
[Roadmap](project/roadmap.md#phase-4--regulatory-submission).

## Not a substitute for standard of care

SmartSepsis-Oph designs do **not** substitute for:

- Bacterial culture.
- Antibiogram.
- Clonal surveillance.
- Any other standard-of-care diagnostic.

## Predicted structures and embeddings

Structures in the [published library](data/library.md) are **predicted**
(ESMFold / ColabFold AF2 / AlphaFold Server AF3) and carry inherent
uncertainty. Use per-variant `struct_mean_plddt` for weighting.

**AlphaFold Server (AF3) terms** apply to the mecA1 / mecA2 PDBs — research
use only, no commercial redistribution.

## Conflict of interest

Authors declare affiliation with **IA para Médicos** (project sponsor). No
financial conflicts.

## Contact

For any questions about scope or regulatory status, write to
**iaparamedicos@gmail.com**.
