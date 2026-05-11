# Manuscript Skeleton — Scoping Review of CRISPR-Based Diagnostics for Ocular Bacterial Infections

**Working title:** CRISPR-Based Diagnostics for Ocular Bacterial Infections: A PRISMA-ScR Scoping Review and Evidence Map for an AI-First Point-of-Care Platform.

**Authors (provisional):** Pacheco Dias JV; Sakuno G; [additional co-authors TBD]
**Corresponding author:** João Victor Pacheco Dias (IA para Médicos / HC-FMUSP).

**Paper 02 of the SmartSepsis-Oph research line.** This skeleton is to be filled after data extraction is complete. No results are pre-stated.

---

## Target journals (in priority order)

1. **Survey of Ophthalmology** — high-impact ophthalmology review venue, accepts scoping reviews; strong fit for evidence-map framing.
2. **BMJ Open Ophthalmology** — open access, methodologically permissive of PRISMA-ScR, good for protocol-paired publication.
3. **Frontiers in Cellular and Infection Microbiology** — fits the molecular-diagnostic angle, especially if AMR markers prominent.
4. **Diagnostics (MDPI)** — fast turnaround, fits if framing leans toward technology landscape and target product profile.

Backup: *Ophthalmology Science*, *Journal of Clinical Microbiology*, *Translational Vision Science & Technology*.

---

## Abstract (structured, 250–300 words)

- Background
- Objective
- Methods (PRISMA-ScR; databases; date range; PCC; two-reviewer screening; OSF registration ID)
- Results (n records identified, n included, primary findings of the evidence map — to be filled)
- Conclusions and implications for an AI-first ophthalmic CRISPR-Dx platform

**Keywords:** CRISPR; Cas12a; SHERLOCK; DETECTR; endophthalmitis; microbial keratitis; point-of-care; isothermal amplification; scoping review; PRISMA-ScR.

---

## 1. Introduction

1.1 Clinical burden of ocular bacterial infections (endophthalmitis, microbial keratitis, ocular surface infection); time-to-pathogen-identification as a determinant of visual outcome.
1.2 Limitations of current diagnostics: culture sensitivity in vitreous and aqueous taps, Gram stain yield in corneal scrapings, turnaround of 16S PCR and mNGS, cost and infrastructure of multiplex panels.
1.3 Emergence of CRISPR-Dx (SHERLOCK 2017, DETECTR 2018) and its theoretical fit for ophthalmic POC use: micro-volume samples, attomolar sensitivity, isothermal pre-amplification, simple readout.
1.4 Rationale for a scoping review: confirm and characterize the evidence gap; map adjacent feasibility evidence; inform an AI-first Cas12a platform (SmartSepsis-Oph).
1.5 Objective and review question (verbatim from protocol).

## 2. Methods

2.1 Protocol and registration (OSF DOI to insert).
2.2 PRISMA-ScR reporting compliance; PRESS-reviewed search strings.
2.3 Eligibility (PCC summary; full criteria in supplement).
2.4 Information sources and search execution dates.
2.5 Selection process (two independent reviewers, kappa on calibration set, third reviewer for conflicts).
2.6 Data charting items (link to protocol Table).
2.7 Synthesis approach (descriptive + tabular + 2D evidence map).
2.8 Deviations from protocol (if any).

## 3. Results

3.1 **PRISMA flow diagram** — to be inserted as Figure 1. Description: records identified (per source), after duplicates removed, screened on title/abstract, sought for retrieval, assessed for eligibility, included; exclusion reasons quantified at full-text stage. Adjacent searches reported in a parallel sub-flow.

3.2 Characteristics of included studies (Table 1).

3.3 Target conditions and pathogens covered (Table 2).

3.4 Technology stack (Cas enzyme × amplification × readout) (Table 3).

3.5 Validation stage distribution (Table 4) and evidence map (Figure 2: pathogen × validation-stage heatmap).

3.6 Adjacent evidence:
- 3.6.1 CRISPR-Dx benchmarks from non-ocular bacterial infection (Adjacent A).
- 3.6.2 Non-CRISPR molecular POC validated on ocular samples (Adjacent B).

### Pre-specified tables

**Table 1 — Included primary studies.** Columns: Study ID; Country; Design; Population / n; Sample type; Cas enzyme; Amplification; Readout; Target pathogen(s) / gene(s); LoD; Sensitivity / Specificity vs reference standard; Time-to-result; Validation stage; Notes.

**Table 2 — Target conditions and pathogens.** Rows: ocular condition (endophthalmitis post-cataract, post-injection, endogenous; microbial keratitis bacterial; bacterial conjunctivitis; orbital cellulitis; other). Columns: pathogens covered by at least one included study (*S. aureus*, CoNS, *S. pneumoniae*, *P. aeruginosa*, *Enterococcus* spp., *Streptococcus* spp., Enterobacterales, *Moraxella*, AMR markers); cell = study IDs.

**Table 3 — Technology stack used by included studies.** Rows: studies. Columns: Cas enzyme; amplification chemistry; readout modality; hardware footprint; AI/ML component (yes/no, description); end-to-end time; consumable cost (if reported).

**Table 4 — Validation stage of included studies.** Rows: studies. Columns: in silico → analytical (synthetic targets) → contrived clinical (spiked ocular matrix) → retrospective clinical samples → prospective clinical → regulatory submission. Color-coded for the figure version.

**Figure 1.** PRISMA-ScR flow diagram (core search + adjacent searches as parallel arms).
**Figure 2.** Evidence map — pathogen × validation stage heatmap; cell color = study count, cell text = platform names.

## 4. Discussion

4.1 **Principal findings.** Summarize size and shape of the evidence base; expected to confirm a substantial gap for ocular bacterial CRISPR-Dx.

4.2 **Gaps identified.** Pathogen × sample-type × validation-stage combinations with zero or near-zero coverage; absence (or presence) of prospective clinical validation in vitreous, aqueous, and corneal scrape matrices; absence (or presence) of integrated AI/ML decision support.

4.3 **Implications for an ophthalmic CRISPR-Dx program.** What the adjacent CRISPR-Dx literature (Adjacent A) suggests is technically achievable for LoD, time-to-result, and AMR multiplexing; what the ophthalmic molecular-POC literature (Adjacent B) suggests is the performance bar to clear and the workflow constraints (sample volume, point of care in OR / clinic / emergency).

4.4 **Where SmartSepsis-Oph fits.** Position the program explicitly and conservatively against the identified gap: an AI-first Cas12a platform targeting the highest-impact ocular bacterial pathogens with pre-specified clinical workflows (post-injection endophthalmitis triage; community microbial keratitis triage in low-resource settings). State explicitly that this is a hypothesis-generating positioning, not an efficacy claim. Avoid overstating: no clinical performance is implied by this review.

4.5 **Comparison with prior reviews.** Narrative comparison to any prior scoping or narrative reviews on (a) CRISPR-Dx broadly, (b) molecular diagnosis of endophthalmitis / keratitis.

4.6 **Strengths of this review.**
- Pre-registered protocol on OSF.
- Two-reviewer independent screening with calibration and third-reviewer adjudication.
- PRESS-peer-reviewed searches across five bibliographic databases, two preprint servers, two trial registries, grey literature, and four society conference archives.
- Pre-specified adjacent triangulation searches to characterize feasibility despite expected low primary yield.
- Multilingual screening (English, Portuguese, Spanish).

4.7 **Limitations of this review.**
- Excluded languages beyond EN/PT/ES.
- CRISPR-Dx evolves rapidly; preprint and conference content may be superseded by the time of publication.
- Heterogeneity of reference standards and validation reporting in primary studies limits cross-study comparison.
- No formal risk-of-bias appraisal (PRISMA-ScR does not require one); feasibility-relevant tags applied instead.
- Conference abstracts may underreport methods and metrics.

4.8 **Future research recommendations.** Specific, scoped, and tied to identified gaps; include a proposed target product profile sketch for an ophthalmic CRISPR-Cas12a platform.

## 5. Conclusions

Concise (1 paragraph). State the evidence-map result and the implication for an AI-first ophthalmic CRISPR-Dx development program. Do not over-claim.

---

## Declarations

- **Funding:** to insert.
- **Conflicts of interest:** to insert per ICMJE; declare any relationship of any author to a CRISPR-Dx or ophthalmic-diagnostic commercial entity.
- **Author contributions:** CRediT taxonomy.
- **Data availability:** all search logs, screening decisions, charting forms, and analysis code archived on the OSF project (DOI to insert). No participant-level data.
- **Ethics:** scoping review of published literature; IRB / CEP not required.
- **AI use disclosure:** disclose any LLM assistance for drafting per target-journal policy.
- **Acknowledgements:** to insert.

## References

To be assembled in target journal style. Reference manager: Zotero with the SmartSepsis-Oph group library.

## Supplementary materials

- S1: Full search strings as executed (verbatim from `scoping_review_search_strategy.md`), with date, hit counts, and platform versions.
- S2: PRISMA-ScR checklist.
- S3: PRESS checklist.
- S4: Calibration screening kappa report.
- S5: Charting form (blank template).
- S6: List of excluded full-text studies with reasons.
- S7: OSF protocol registration record.
