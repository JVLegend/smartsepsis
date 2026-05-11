# Search Strategy — CRISPR-Dx for Ocular Bacterial Infections (PRISMA-ScR)

Companion to `scoping_review_protocol.md`. All strings are executable as written; date and result counts to be logged at execution time. Search window: **2017-01-01 to date of execution**. No language restriction at the database query level (language filter applied at screening). PRESS 2015 peer review of these strings to be completed before execution.

---

## A. Core search — CRISPR-Dx applied to ocular bacterial infection

Three concept blocks combined with AND:
- **Block 1 (CRISPR-Dx)**
- **Block 2 (Diagnostic / point-of-care / isothermal)**
- **Block 3 (Eye / ocular / ophthalmic)**

### A.1 PubMed / MEDLINE

```
(
  "CRISPR-Cas Systems"[MeSH] OR "Clustered Regularly Interspaced Short Palindromic Repeats"[MeSH]
  OR CRISPR[tiab] OR "Cas9"[tiab] OR "Cas12"[tiab] OR "Cas12a"[tiab] OR "Cas12b"[tiab]
  OR "Cas13"[tiab] OR "Cas13a"[tiab] OR "Cas13b"[tiab] OR "Cas14"[tiab]
  OR SHERLOCK[tiab] OR DETECTR[tiab] OR HOLMES[tiab] OR CDetection[tiab] OR FELUDA[tiab]
  OR "collateral cleavage"[tiab] OR "trans-cleavage"[tiab] OR "guide RNA"[tiab] OR crRNA[tiab]
)
AND
(
  "Molecular Diagnostic Techniques"[MeSH] OR "Point-of-Care Testing"[MeSH] OR "Nucleic Acid Amplification Techniques"[MeSH]
  OR diagnos*[tiab] OR detection[tiab] OR "point of care"[tiab] OR "point-of-care"[tiab]
  OR "nucleic acid test"[tiab] OR NAAT[tiab]
  OR "isothermal amplification"[tiab] OR "recombinase polymerase amplification"[tiab] OR RPA[tiab]
  OR "loop-mediated isothermal amplification"[tiab] OR LAMP[tiab] OR "RT-LAMP"[tiab]
  OR "lateral flow"[tiab]
)
AND
(
  "Eye"[MeSH] OR "Eye Infections, Bacterial"[MeSH] OR "Endophthalmitis"[MeSH] OR "Keratitis"[MeSH]
  OR "Conjunctivitis, Bacterial"[MeSH] OR "Uveitis"[MeSH] OR "Cornea"[MeSH]
  OR "Aqueous Humor"[MeSH] OR "Vitreous Body"[MeSH] OR "Tears"[MeSH]
  OR eye[tiab] OR ocular[tiab] OR ophthalm*[tiab] OR endophthalmitis[tiab]
  OR keratitis[tiab] OR conjunctivitis[tiab] OR uveitis[tiab]
  OR vitreous[tiab] OR cornea*[tiab] OR "aqueous humor"[tiab] OR "ocular surface"[tiab]
  OR tear[tiab] OR tears[tiab] OR "contact lens"[tiab]
)
AND ("2017/01/01"[Date - Publication] : "3000"[Date - Publication])
```

Filter: Humans NOT applied at query level (analytical and animal-model studies are eligible). No language filter.

### A.2 Embase (Elsevier, emtree.com syntax)

```
(
  'crispr cas system'/exp OR crispr:ti,ab,kw OR 'cas9':ti,ab,kw OR 'cas12':ti,ab,kw
  OR 'cas12a':ti,ab,kw OR 'cas13':ti,ab,kw OR 'cas14':ti,ab,kw
  OR sherlock:ti,ab,kw OR detectr:ti,ab,kw OR holmes:ti,ab,kw OR feluda:ti,ab,kw
  OR 'collateral cleavage':ti,ab,kw OR 'trans cleavage':ti,ab,kw
  OR 'guide rna':ti,ab,kw OR crrna:ti,ab,kw
)
AND
(
  'molecular diagnosis'/exp OR 'point of care testing'/exp OR 'nucleic acid amplification'/exp
  OR diagnos*:ti,ab,kw OR detection:ti,ab,kw OR 'point of care':ti,ab,kw
  OR 'isothermal amplification':ti,ab,kw OR rpa:ti,ab,kw OR lamp:ti,ab,kw
  OR 'rt lamp':ti,ab,kw OR 'lateral flow':ti,ab,kw OR 'nucleic acid test':ti,ab,kw
)
AND
(
  'eye'/exp OR 'eye infection'/exp OR 'endophthalmitis'/exp OR 'keratitis'/exp
  OR 'bacterial conjunctivitis'/exp OR 'uveitis'/exp OR 'cornea'/exp
  OR 'aqueous humor'/exp OR 'vitreous body'/exp OR 'tear'/exp
  OR eye:ti,ab,kw OR ocular:ti,ab,kw OR ophthalm*:ti,ab,kw
  OR endophthalmitis:ti,ab,kw OR keratitis:ti,ab,kw OR conjunctivitis:ti,ab,kw
  OR uveitis:ti,ab,kw OR vitreous:ti,ab,kw OR cornea*:ti,ab,kw
  OR 'aqueous humor':ti,ab,kw OR 'ocular surface':ti,ab,kw
  OR tear:ti,ab,kw OR tears:ti,ab,kw OR 'contact lens':ti,ab,kw
)
AND [2017-2026]/py
```

### A.3 Scopus

```
TITLE-ABS-KEY(
  ( CRISPR OR "Cas9" OR "Cas12" OR "Cas12a" OR "Cas12b" OR "Cas13" OR "Cas13a" OR "Cas13b" OR "Cas14"
    OR SHERLOCK OR DETECTR OR HOLMES OR CDetection OR FELUDA
    OR "collateral cleavage" OR "trans-cleavage" OR "guide RNA" OR crRNA )
)
AND TITLE-ABS-KEY(
  ( diagnos* OR detection OR "point of care" OR "point-of-care"
    OR "nucleic acid test" OR NAAT
    OR "isothermal amplification" OR "recombinase polymerase amplification" OR RPA
    OR "loop-mediated isothermal amplification" OR LAMP OR "RT-LAMP"
    OR "lateral flow" )
)
AND TITLE-ABS-KEY(
  ( eye OR ocular OR ophthalm* OR endophthalmitis OR keratitis OR conjunctivitis OR uveitis
    OR vitreous OR cornea* OR "aqueous humor" OR "ocular surface" OR tear OR tears OR "contact lens" )
)
AND PUBYEAR > 2016
```

### A.4 Web of Science Core Collection

```
TS=(
  ( CRISPR OR "Cas9" OR "Cas12" OR "Cas12a" OR "Cas13" OR "Cas13a" OR "Cas14"
    OR SHERLOCK OR DETECTR OR HOLMES OR FELUDA
    OR "collateral cleavage" OR "trans-cleavage" OR "guide RNA" OR crRNA )
)
AND TS=(
  ( diagnos* OR detection OR "point of care" OR "point-of-care"
    OR "nucleic acid test" OR NAAT
    OR "isothermal amplification" OR "recombinase polymerase amplification" OR RPA
    OR "loop-mediated isothermal amplification" OR LAMP OR "RT-LAMP" OR "lateral flow" )
)
AND TS=(
  ( eye OR ocular OR ophthalm* OR endophthalmitis OR keratitis OR conjunctivitis OR uveitis
    OR vitreous OR cornea* OR "aqueous humor" OR "ocular surface" OR tear OR tears OR "contact lens" )
)
AND PY=(2017-2026)
```

### A.5 Cochrane CENTRAL

Same Block 1 × Block 2 × Block 3 structure, using CENTRAL's simplified search syntax with MeSH descriptors where available. Limit: 2017–present.

### A.6 Preprint servers (bioRxiv, medRxiv)

Use the Europe PMC preprint interface (https://europepmc.org/) with the PubMed string A.1 filtered to `SRC:PPR`, plus native bioRxiv/medRxiv search using the simplified string:
```
(CRISPR OR Cas12 OR Cas13 OR SHERLOCK OR DETECTR) AND (diagnostic OR detection OR point-of-care OR isothermal OR LAMP OR RPA) AND (eye OR ocular OR ophthalm* OR endophthalmitis OR keratitis OR conjunctivitis OR vitreous OR cornea)
```

### A.7 Trial registries

- **ClinicalTrials.gov advanced search:** Condition/disease = (endophthalmitis OR keratitis OR conjunctivitis OR "ocular infection" OR "eye infection"); Other terms = (CRISPR OR Cas12 OR Cas13 OR SHERLOCK OR DETECTR OR isothermal OR LAMP OR RPA); Start date 2017-01-01.
- **WHO ICTRP:** same terms in basic search; export and screen.

### A.8 Google Scholar (grey literature)

Run as two queries; first 200 results of each, sorted by relevance, screened in order:

1. `(CRISPR OR Cas12 OR Cas13 OR SHERLOCK OR DETECTR) (endophthalmitis OR keratitis OR "ocular infection")`
2. `("point of care" OR isothermal OR LAMP OR RPA) CRISPR (eye OR ocular OR vitreous OR cornea)`

Records exported via Publish-or-Perish for reproducibility; date and total hit count logged.

### A.9 Conference abstracts

Hand-searched in the official abstract archives for 2017–present:
- **ARVO Annual Meeting** (IOVS abstract supplements) — search terms: CRISPR, Cas12, Cas13, SHERLOCK, DETECTR, isothermal, LAMP, RPA.
- **AAO Subspecialty Day and Annual Meeting** abstracts.
- **ESCRS, APAO, WOC** abstract archives.

---

## B. Adjacent triangulation searches (pre-specified, tagged separately)

These searches are **not** part of the core PCC map but are pre-specified to triangulate feasibility and benchmarks. Reported in a dedicated subsection of the manuscript.

### B.1 Adjacent search A — CRISPR-Dx for bacterial infection in non-ocular body sites

Used to establish performance benchmarks of CRISPR-Dx against bacterial pathogens generally. Run on PubMed and Scopus only.

```
(CRISPR OR Cas12 OR Cas13 OR SHERLOCK OR DETECTR OR HOLMES)
AND (diagnos* OR detection OR "point of care" OR isothermal OR RPA OR LAMP)
AND (sepsis OR bloodstream OR bacteremia OR pneumonia OR respiratory OR urinary tract OR meningitis OR cerebrospinal OR tuberculosis OR "antimicrobial resistance" OR mecA OR vanA OR KPC OR NDM)
AND PUBYEAR > 2016
```

### B.2 Adjacent search B — molecular POC and isothermal amplification on ocular samples (non-CRISPR)

Used to map the current ophthalmic molecular-POC landscape that a CRISPR-Dx product would have to outperform. Run on PubMed and Scopus.

```
( "point of care" OR "point-of-care" OR isothermal OR LAMP OR "RT-LAMP" OR RPA
  OR "multiplex PCR" OR "metagenomic sequencing" OR mNGS OR "nanopore" )
AND ( endophthalmitis OR keratitis OR conjunctivitis OR uveitis OR "ocular infection"
      OR "eye infection" OR vitreous OR aqueous OR cornea* OR "ocular surface" )
AND NOT (CRISPR OR Cas12 OR Cas13)
AND PUBYEAR > 2014
```

(Adjacent B uses a wider 2015+ window because non-CRISPR molecular POC predates 2017.)

---

## C. Expected yield and rationale

The core search (Section A) is expected to return a **small number of primary records** (plausibly fewer than ~30 unique studies after deduplication, with substantially fewer meeting full eligibility) because:

1. CRISPR-Dx as a field begins effectively in 2017.
2. Ocular bacterial infection is a narrow clinical niche with limited sample volume and limited industry incentive.
3. Most CRISPR-Dx development has prioritized respiratory viruses (SARS-CoV-2), bloodstream, and tropical pathogens.

This expected low yield is **the principal finding** of the scoping review and justifies the adjacent searches (Section B), which will be reported with clear methodological separation. Total volume across all searches is expected to be tractable for two-reviewer screening within the protocol timeline.

---

## D. Record management

- Export format: RIS / NBIB / CSV depending on source.
- Deduplication: Zotero + manual reconciliation on DOI, then title + first-author + year fuzzy match.
- Screening platform: Rayyan (preferred) or Covidence.
- All raw export files, dedup logs, screening decisions, and conflict-resolution notes archived to the OSF project before manuscript submission.
