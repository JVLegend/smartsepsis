# PRISMA-ScR Protocol: CRISPR-Based Diagnostics for Ocular Bacterial Infections

**Paper 02 of the SmartSepsis-Oph research line**

- **Lead reviewer:** João Victor Pacheco Dias (IA para Médicos / HC-FMUSP doctoral candidate)
- **Clinical advisor / second reviewer:** Dr. Gustavo Sakuno (Mass Eye and Ear / Harvard Medical School)
- **Third reviewer (conflict resolution):** to be designated prior to screening
- **Protocol version:** 0.1 (draft, pre-registration)
- **Date:** 2026-05-11

---

## 1. Title

A scoping review of CRISPR-based nucleic acid diagnostics for ocular bacterial infections and adjacent point-of-care molecular technologies relevant to an AI-first Cas12a platform (PRISMA-ScR).

## 2. Background and Rationale

Bacterial infections of the eye — endophthalmitis, microbial keratitis, and severe ocular surface infections — are sight-threatening emergencies in which time-to-pathogen identification is a primary determinant of visual outcome. Standard-of-care microbiology (culture, Gram stain, broad-range 16S PCR on aqueous or vitreous taps) is slow, often low-yield, and rarely available at the bedside or in low-resource settings.

CRISPR-based nucleic acid detection platforms (SHERLOCK using Cas13, DETECTR using Cas12a, and derivatives) emerged from 2017 onward and have demonstrated attomolar sensitivity, single-nucleotide specificity, and compatibility with isothermal pre-amplification (RPA, LAMP, RT-LAMP) and lateral-flow or fluorescence readout. These properties make CRISPR-Dx an attractive substrate for point-of-care (POC) molecular diagnosis in ophthalmology, where sample volumes are micro-liter scale and turnaround time is clinically decisive.

To our knowledge, very little primary literature describes CRISPR-Dx applied directly to ocular bacterial pathogens or ocular biological matrices. A rigorous scoping review is therefore needed to (i) confirm and characterize this evidence gap, (ii) map adjacent evidence (CRISPR-Dx for bacterial infection in other anatomical sites, molecular POC platforms validated on ocular samples, isothermal amplification applied to eye fluids), and (iii) inform the design of an AI-first Cas12a platform for ophthalmic infection (SmartSepsis-Oph, Paper 03+).

## 3. Objective and Review Question

**Primary objective:** Map the published and grey literature on CRISPR-based diagnostic development for ocular bacterial infections, and characterize adjacent evidence relevant to translation.

**Review question:** What is the current state of CRISPR-based diagnostic development for ocular bacterial infections (endophthalmitis, microbial keratitis, ocular surface infections), and what adjacent evidence informs feasibility of an AI-first point-of-care CRISPR-Cas12a platform for these conditions?

**Sub-questions:**
1. How many primary studies report CRISPR-Dx evaluated on ocular pathogens or ocular biological samples? What is their validation stage (in silico, analytical, contrived clinical, prospective clinical)?
2. Which Cas enzymes, amplification chemistries, target genes, readouts, and limits of detection have been reported?
3. What adjacent molecular POC platforms (multiplex PCR panels, isothermal amplification, metagenomic NGS) have been validated on eye samples, and what do their performance benchmarks suggest for a CRISPR-Dx target product profile?
4. What gaps exist that an AI-first CRISPR-Cas12a ophthalmic platform could plausibly address?

## 4. Eligibility Criteria (PCC Framework)

**Population (P):**
- Human patients with, or at risk of, ocular bacterial infection (endophthalmitis post-cataract / post-injection / endogenous, microbial keratitis, infectious conjunctivitis, blepharitis, orbital cellulitis, ocular surface infection).
- Ocular biological samples: aqueous humor, vitreous humor, corneal scrapings or buttons, conjunctival swabs, tear film, eyelid swabs, contact lens cases / solutions used as proxies.
- Animal models of ocular bacterial infection accepted as adjacent evidence (tagged separately).

**Concept (C):**
- CRISPR-based nucleic acid diagnostics: Cas9, Cas12 (a/b/f), Cas13 (a/b/d), Cas14, and derivatives.
- Named platforms: SHERLOCK, DETECTR, HOLMES, CDetection, STOPCovid, AIOD-CRISPR, FELUDA, and successors.
- Isothermal amplification coupled to CRISPR readout: RPA, LAMP, RT-LAMP, NASBA, RCA + CRISPR.
- Lateral-flow, fluorescence, electrochemical, or smartphone-based CRISPR readouts.

**Context (C):**
- Clinical, pre-clinical, analytical, or feasibility studies.
- Any country, any healthcare setting (tertiary, community, low-resource, field).
- Languages: English, Portuguese, Spanish.
- Time window: 2017-01-01 to date of search execution (CRISPR-Dx field effectively begins with Gootenberg et al. 2017).

**Study designs included:**
- Primary research: analytical validation studies, diagnostic accuracy studies, contrived-sample studies, prospective and retrospective clinical evaluations, animal model validation, device feasibility / engineering reports.
- Conference abstracts and preprints (tagged as grey literature).
- Registered trial protocols (ClinicalTrials.gov, WHO ICTRP).

**Study designs excluded:**
- Reviews, editorials, commentaries, letters without primary data (retained only for citation chasing).
- CRISPR studies focused exclusively on genome editing or therapeutic use without diagnostic intent.
- Diagnostic studies of non-bacterial ocular pathology where no bacterial target is reported (viral-only, fungal-only, and parasitic-only studies are excluded from the core map but tagged in the adjacent search if they use CRISPR on ocular samples).

## 5. Information Sources

**Bibliographic databases:**
- PubMed / MEDLINE
- Embase (Elsevier)
- Scopus
- Web of Science Core Collection
- Cochrane CENTRAL

**Preprint servers:**
- bioRxiv
- medRxiv

**Trial registries:**
- ClinicalTrials.gov
- WHO ICTRP

**Grey literature:**
- Google Scholar — first 200 results per concept block, screened in order.
- Conference abstracts: ARVO (Annual Meeting), AAO (American Academy of Ophthalmology), ESCRS (European Society of Cataract & Refractive Surgeons), APAO, WOC.
- Society and consortium reports relevant to molecular ophthalmic diagnostics.

**Citation chasing:**
- Forward and backward citation tracking on every included primary study and on any narrative review identified during screening.

## 6. Search Strategy (Concept)

Three Boolean concept blocks combined with AND:

1. **CRISPR-Dx concept:** CRISPR, Cas9, Cas12, Cas13, Cas14, SHERLOCK, DETECTR, HOLMES, collateral cleavage, trans-cleavage, guide RNA, crRNA.
2. **Diagnostic / POC concept:** diagnostic, detection, point-of-care, nucleic acid test, NAAT, isothermal amplification, RPA, LAMP, RT-LAMP, lateral flow.
3. **Ocular / ophthalmic concept:** eye, ocular, ophthalmic, endophthalmitis, keratitis, conjunctivitis, uveitis, vitreous, aqueous humor, cornea, ocular surface, tear, contact lens.

Full executable strings, per database, are provided in `scoping_review_search_strategy.md`. Given expected low primary yield, two **adjacent triangulation searches** are pre-specified:

- **Adjacent A:** CRISPR-Dx concept AND Diagnostic/POC concept AND (sepsis OR bloodstream OR respiratory OR urinary OR meningitis OR cerebrospinal) — to map CRISPR-Dx performance benchmarks in other body compartments.
- **Adjacent B:** Diagnostic/POC concept (excluding CRISPR) AND Ocular/ophthalmic concept AND (PCR OR multiplex OR metagenomic OR isothermal OR LAMP OR mNGS) — to map the molecular POC landscape in ophthalmology that a CRISPR-Dx product would have to outperform.

## 7. Study Selection Process

1. Records from all sources exported to a Zotero group library, deduplicated, then loaded into Rayyan (or Covidence) for blinded screening.
2. **Title / abstract screening:** two reviewers (JVPD, GS) independently screen each record against eligibility criteria. Inter-rater agreement reported as Cohen's kappa on a calibration set of the first 100 records before full screening proceeds.
3. **Full-text screening:** same two reviewers independently retrieve and assess full texts. Reasons for exclusion recorded.
4. **Conflict resolution:** disagreements resolved by discussion; unresolved conflicts adjudicated by a pre-designated third reviewer.
5. PRISMA-ScR flow diagram populated and reported.

## 8. Data Charting / Extraction Template

A piloted charting form (Google Sheets / REDCap) will capture, per included study:

| Field | Description |
|---|---|
| Study ID | First author, year |
| Country | Affiliation of corresponding author / site of clinical samples |
| Design | Analytical / contrived / retrospective clinical / prospective clinical / animal / engineering |
| Population | Condition, n, age, setting |
| Sample type | Aqueous, vitreous, corneal scrape, conjunctival swab, tear, contact lens, other |
| Target pathogen(s) | Species and resistance markers (mecA, vanA, etc.) |
| Target gene(s) | 16S rRNA, species-specific, AMR markers |
| Cas enzyme | Cas9 / Cas12a / Cas12b / Cas13a / Cas13b / Cas14 / other |
| Amplification | None / RPA / LAMP / RT-LAMP / PCR / NASBA |
| Readout | Fluorescence / lateral flow / electrochemical / smartphone / other |
| Limit of detection | Reported LoD with units |
| Analytical sensitivity / specificity | If reported |
| Clinical sensitivity / specificity | If reported, with reference standard |
| Time-to-result | Minutes, end-to-end |
| Reference standard | Culture / 16S PCR / mNGS / clinical adjudication |
| Validation stage | TRL-style: in silico, analytical, contrived clinical, prospective clinical, regulatory |
| Funding / COI | As declared |
| Notes | AI / ML components, hardware footprint, deployment context |

Charting will be piloted on the first 5 included studies independently by both reviewers; the form will be revised iteratively per PRISMA-ScR guidance.

## 9. Synthesis Approach

- **Descriptive synthesis** of study characteristics, technology stack, target pathogens, and validation maturity.
- **Tabular summary** (Tables 1–4 in `scoping_review_skeleton.md`).
- **Evidence map** as a 2D heatmap: pathogen / condition (rows) by validation stage (columns), with cell color encoding number of studies and cell content listing platform(s).
- **Narrative gap analysis** anchored to the PCC framework, with explicit identification of unstudied pathogen × sample-type × platform combinations.
- No meta-analysis (scoping review; quality appraisal not required by PRISMA-ScR but a light feasibility-relevant signal — e.g., whether prospective clinical samples were used — will be tagged).

## 10. Registration

This protocol will be registered on the Open Science Framework (OSF Registries) prior to the start of data extraction. The OSF registration DOI will be inserted here once issued. PROSPERO does not accept scoping reviews; OSF is the appropriate registry.

## 11. Ethics

Scoping reviews of published and grey literature do not involve human subjects research and do not require IRB / CEP review. No identifiable patient data will be handled. If individual-patient data appear in any included abstract or preprint, they will not be re-extracted at the participant level.

## 12. Timeline (rough, in months from kickoff)

| Month | Milestone |
|---|---|
| 0 | Protocol finalized, OSF registration submitted |
| 1 | Search strings peer-reviewed (PRESS checklist), searches executed across all databases and grey sources |
| 1–2 | Deduplication, calibration screening (kappa), title/abstract screening |
| 2–3 | Full-text screening and conflict resolution |
| 3–4 | Data charting (piloted on 5, then completed) |
| 4–5 | Synthesis, evidence map construction, table assembly |
| 5–6 | Manuscript drafting (skeleton in `scoping_review_skeleton.md`) |
| 6 | Internal review with clinical advisor, submission to target journal |

## 13. Reporting

Reporting will follow the PRISMA-ScR checklist (Tricco et al., Ann Intern Med 2018). Search strategies will additionally be peer-reviewed against the PRESS 2015 checklist before execution.

## 14. Amendments

Any protocol amendments after OSF registration will be logged with date, rationale, and reviewer initials in an amendments appendix of the final manuscript.
