# SmartSepsis-Oph

**AI-first molecular diagnostics design for ophthalmic infections.**

!!! warning "Phase 0 — computational design only"
    SmartSepsis-Oph is a research program in **computational design phase**.
    All results in this repository are *in silico*. **No claims of clinical
    performance are made.** Experimental validation is pending and is what we
    are actively seeking collaborators for. This is not a medical device and
    has not been approved by ANVISA, FDA, or any regulatory body.

Bacterial endophthalmitis can blind an eye in 24–48 hours. Cultures take
48–72. Existing molecular diagnostics are built for sample volumes orders of
magnitude larger than what ophthalmology yields (sub-microliter vitreous tap,
corneal scrape, ocular swab).

SmartSepsis-Oph attacks this gap with an AI-first design pipeline:

- **Guide-RNA library design** (CRISPR-Cas12a) for 12 resistance-gene families
  relevant to ocular infection, prioritized for the Brazilian epidemiological
  context.
- **In silico specificity modeling** against public reference repositories.
- **Isothermal RPA primer design** compatible with thermocycler-free,
  point-of-care workflows.
- **Functional impact scoring** of variants (inspired by Evo 2 / EVEE).
- **Paper-strip assay architecture** under design.

First declared target: ***Staphylococcus aureus* / mecA** — the most common
and best-validated combination in the CRISPR-Dx literature.

## Quick links

<div class="grid cards" markdown>

-   :material-microscope: **[Methods](methods/pipeline-overview.md)**

    The seven-step design pipeline, from sequence fetch to multiplex panel.

-   :material-database: **[Data](data/library.md)**

    The published multimodal AMR library on Hugging Face.

-   :material-map: **[Roadmap](project/roadmap.md)**

    Phased plan from in silico design to regulatory submission.

-   :material-handshake: **[Partners](project/partners.md)**

    What we are seeking — wet-lab, advisors, clinical, funding.

-   :material-source-pull: **[Contributing](contributing.md)**

    How to engage, propose changes, or join the team.

-   :material-format-quote-close: **[Citing](citing.md)**

    Citation metadata for the software and the dataset.

</div>

## Status at a glance

| Phase | Scope | Status |
|---|---|---|
| **0 — Computational design** | Guide and primer libraries, in silico specificity, scoring, structural analysis. | **Current** |
| **1 — Wet-lab feasibility** | Single-organism, single-gene validation (S. aureus / mecA). | Seeking wet-lab partner |
| **2 — Panel expansion** | Multiplex of priority resistance markers. | Planned |
| **3 — Clinical validation** | Prospective sampling under IRB. | Planned |
| **4 — Regulatory submission** | IVD (RDC 830/2023) pathway and equivalent. | Planned |

See the full [roadmap](project/roadmap.md).

## Contact

- Email: **iaparamedicos@gmail.com**
- Website: <https://smartsepsis.vercel.app>
- Issues: <https://github.com/JVLegend/smartsepsis/issues>

Co-authorship on the forthcoming framework preprint is offered to early
collaborators. See [Partners](project/partners.md).
