# Roadmap

!!! info "Honest framing"
    Status badges below reflect what has actually been delivered. We do not
    claim progress in phases we have not yet entered.

## Phase 0 — Computational design  ![status](https://img.shields.io/badge/status-current-orange)

**Scope.** Guide and primer libraries, in silico specificity, functional
scoring, structural analysis.

**Deliverables:**

- [x] CRISPR-Cas12a guide-RNA design pipeline ([guide design](../methods/guide-design.md)).
- [x] Covariance probe re-scoring with 18 biophysical features ([scoring](../methods/scoring.md)).
- [x] Isothermal RPA primer design.
- [x] In silico specificity check via BLAST ([specificity](../methods/specificity.md)).
- [x] CARD database enrichment.
- [x] Published multimodal AMR dataset on Hugging Face ([library](../data/library.md)).
- [ ] Framework preprint.
- [ ] Migration of pipeline scripts under a packaged `src/` layout.

## Phase 1 — Wet-lab feasibility  ![status](https://img.shields.io/badge/status-seeking%20partner-red)

**Scope.** Single-organism, single-gene validation (***S. aureus* / mecA**) in
cultured isolates.

**Deliverables (planned):**

- [ ] Cas12a collateral-cleavage validation on synthetic targets.
- [ ] RPA amplification check on bacterial genomic DNA.
- [ ] Specificity stress-test against close non-target genomes.
- [ ] Limit-of-detection envelope (preliminary).

## Phase 2 — Panel expansion  ![status](https://img.shields.io/badge/status-planned-lightgrey)

**Scope.** Multiplex of priority resistance markers and real-amplicon validation.

- [ ] Multiplex panel optimization ([multiplex](../methods/multiplex.md)).
- [ ] Cross-reactivity mapping inside the panel.
- [ ] Paper-strip / lateral flow prototype.

## Phase 3 — Clinical validation  ![status](https://img.shields.io/badge/status-planned-lightgrey)

**Scope.** Prospective sampling under IRB; head-to-head against culture and PCR.

- [ ] IRB protocol with clinical partner.
- [ ] Sample-handling SOPs for sub-microliter ocular fluids.
- [ ] Prospective performance evaluation.

## Phase 4 — Regulatory submission  ![status](https://img.shields.io/badge/status-planned-lightgrey)

**Scope.** IVD pathway (ANVISA RDC 830/2023) and equivalent.

- [ ] Quality management system.
- [ ] Analytical and clinical performance dossier.
- [ ] Regulatory submission.

---

## How to move us forward

If you can help bridge Phase 0 → Phase 1, see [Partners](partners.md) or write
to **iaparamedicos@gmail.com**.
