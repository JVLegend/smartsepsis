# Specificity Stress-Test Operator Manual

This document explains how to run an in silico off-target / cross-reactivity
stress-test on the CRISPR-Cas12a guide RNAs designed by the SmartSepsis-Oph
pipeline using `scripts/stress_test_specificity.py`.

> **In silico disclaimer.** Every result produced by this workflow is
> computational. A guide that "passes" all three BLAST screens is a
> *candidate* for wet-lab validation, not a clinically specific reagent.
> No analytical or clinical specificity claim can be made from BLAST
> alone; orthogonal wet-lab cross-reactivity panels are still required.

---

## 1. Prerequisites

- NCBI BLAST+ command-line tools (`blastn`, `makeblastdb`) in `$PATH`.
- Python 3.10+ with `biopython` (already in `requirements.txt`).
- `datasets` CLI from NCBI (optional, convenient for genome download), or
  plain `wget` + `gunzip`. `efetch` (Entrez Direct) is used for small
  targeted RefSeq pulls.
- Roughly 5-10 GB of free disk for the human genome BLAST DB; the
  microbiome and fungi DBs are small (< 1 GB combined).

All commands below assume you are working from the repo root and storing
references under `refs/`:

```bash
mkdir -p refs/{human,microbiome,fungi_acanthamoeba} blast_dbs results
```

---

## 2. Reference data to download

### (a) Human genome — GRCh38

Source: NCBI RefSeq, assembly `GCF_000001405.40` (GRCh38.p14).

```bash
# Using NCBI datasets CLI (recommended)
datasets download genome accession GCF_000001405.40 \
    --include genome --filename refs/human/grch38.zip
unzip -j refs/human/grch38.zip "ncbi_dataset/data/GCF_000001405.40/*.fna" \
    -d refs/human/

# Or plain FTP
wget -P refs/human/ \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip refs/human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

### (b) Ocular surface microbiome

Two complementary references are recommended; combine into a single FASTA.

1. **SILVA 138.1 SSU NR99** (16S rRNA), as a broad bacterial backbone,
   subsetted to genera reported in healthy ocular surface studies
   (e.g. *Corynebacterium*, *Staphylococcus epidermidis*, *Cutibacterium
   acnes*, *Streptococcus*, *Moraxella*, *Pseudomonas*, *Bacillus*).

   ```bash
   wget -P refs/microbiome/ \
     https://data.arb-silva.de/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
   gunzip refs/microbiome/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
   ```

2. **Published ocular-surface microbiome reads / contigs.** Recommended
   BioProjects (cite the underlying studies in your preprint methods):

   - PRJNA646695 - Ocular surface 16S, healthy vs. dry-eye (Doan et al.
     style cohorts).
   - PRJNA428163 - Ocular microbiome metagenomics (Doan et al.,
     *Translational Vision Science & Technology*).
   - PRJEB29302 - European ocular-surface 16S cohort.

   Pull representative assemblies / reads:

   ```bash
   datasets download genome taxon "Corynebacterium" --reference \
       --filename refs/microbiome/coryne.zip
   datasets download genome taxon "Staphylococcus epidermidis" --reference \
       --filename refs/microbiome/sepi.zip
   datasets download genome taxon "Cutibacterium acnes" --reference \
       --filename refs/microbiome/cacnes.zip
   datasets download genome taxon "Moraxella" --reference \
       --filename refs/microbiome/moraxella.zip
   # unzip and concatenate all .fna files
   for z in refs/microbiome/*.zip; do unzip -o "$z" -d refs/microbiome/_unz; done
   find refs/microbiome/_unz -name "*.fna" -exec cat {} + \
       > refs/microbiome/ocular_microbiome.fasta
   cat refs/microbiome/SILVA_138.1_SSURef_NR99_tax_silva.fasta \
       >> refs/microbiome/ocular_microbiome.fasta
   ```

   Suggested citations to include in the preprint methods:
   - Doan T. et al., "Paucibacterial Microbiome and Resident DNA Virome of
     the Healthy Conjunctiva." *IOVS* 2016.
   - Ozkan J. et al., "Temporal Stability and Composition of the Ocular
     Surface Microbiome." *Sci Rep* 2017.
   - Petrillo F. et al., "Current Evidence on the Ocular Surface
     Microbiota and Related Diseases." *Microorganisms* 2020.

### (c) Acanthamoeba + ocular fungal negative controls

```bash
# Acanthamoeba castellanii (RefSeq)
datasets download genome accession GCF_000313135.1 \
    --include genome --filename refs/fungi_acanthamoeba/acastellanii.zip

# Aspergillus fumigatus Af293
datasets download genome accession GCF_000002655.1 \
    --include genome --filename refs/fungi_acanthamoeba/afumigatus.zip

# Candida albicans SC5314
datasets download genome accession GCF_000182965.3 \
    --include genome --filename refs/fungi_acanthamoeba/calbicans.zip

# Fusarium solani species complex
datasets download genome accession GCF_020744495.1 \
    --include genome --filename refs/fungi_acanthamoeba/fsolani.zip

for z in refs/fungi_acanthamoeba/*.zip; do
    unzip -o "$z" -d refs/fungi_acanthamoeba/_unz
done
find refs/fungi_acanthamoeba/_unz -name "*.fna" -exec cat {} + \
    > refs/fungi_acanthamoeba/ocular_fungi_acanthamoeba.fasta
```

Or, for any single accession, with Entrez Direct:

```bash
efetch -db nuccore -id NC_007194.1 -format fasta \
    > refs/fungi_acanthamoeba/afumigatus_chr1.fasta
```

---

## 3. Build the three BLAST databases

```bash
makeblastdb -in refs/human/GCF_000001405.40_GRCh38.p14_genomic.fna \
    -dbtype nucl -parse_seqids -title "human_genome" \
    -out blast_dbs/human_genome

makeblastdb -in refs/microbiome/ocular_microbiome.fasta \
    -dbtype nucl -parse_seqids -title "ocular_microbiome" \
    -out blast_dbs/ocular_microbiome

makeblastdb -in refs/fungi_acanthamoeba/ocular_fungi_acanthamoeba.fasta \
    -dbtype nucl -parse_seqids -title "ocular_fungi_acanthamoeba" \
    -out blast_dbs/ocular_fungi_acanthamoeba
```

Once built, the stress test does NOT require network access.

---

## 4. Prepare the guide input

The script accepts CSV, TSV, or parquet with these columns:

| column           | description                                |
|------------------|--------------------------------------------|
| `gene_family`    | e.g. `blaKPC`, `mecA`                      |
| `variant`        | e.g. `KPC-2`                               |
| `spacer_sequence`| 20-25 nt protospacer (no PAM, ACGT only)   |
| `strand`         | `+` or `-`                                 |

The pipeline's `*_cas12a_guides.tsv` files (see `design_guides.py`) can
be converted easily by renaming `spacer_seq` to `spacer_sequence` and
adding the `gene_family` / `variant` columns. The HuggingFace dataset
panel parquet (`data/hf_dataset/panel.parquet`) already exposes these
columns and can be used directly.

---

## 5. Run the stress test

```bash
python scripts/stress_test_specificity.py \
    --input data/hf_dataset/panel.parquet \
    --db-human   blast_dbs/human_genome \
    --db-microb  blast_dbs/ocular_microbiome \
    --db-fungi   blast_dbs/ocular_fungi_acanthamoeba \
    --outdir     results/ \
    --threads    8
```

Optional flags:

- `--max-hits-18nt 0`  - allowed hits with >=18 nt complementarity (default 0).
- `--max-identity-pct 85` - best-hit %identity must stay below this.
- `--no-seed-strict`   - stop treating seed (pos 1-8) hits as automatic fails.

Outputs:

- `results/stress_test_<UTC_timestamp>.tsv` - one row per (guide, DB, hit).
- `results/stress_test_summary.tsv` - per-guide pass/fail per DB plus
  `overall_pass`.

---

## 6. Interpreting the results

A guide is reported as `pass=False` for a database when any of:

- it has >= 1 hit of >= 18 nt complementarity (default), or
- best-hit %identity >= 85 (default), or
- a hit engages the Cas12a seed region (spacer positions 1-8) AND
  `--seed-strict` is on (default).

What to do per failure mode:

| Failed against            | Likely concern                                  | Action                                                                                  |
|---------------------------|-------------------------------------------------|-----------------------------------------------------------------------------------------|
| human_genome              | Host-DNA cleavage / background                  | Redesign: shift PAM, prefer guides absent in GRCh38 BLAST; flag in panel summary.       |
| ocular_microbiome         | Commensal cross-reactivity                      | Inspect hit taxonomy; if a conserved commensal hits in the seed, drop the guide.        |
| ocular_fungi_acanthamoeba | Negative-control reagent has unexpected match   | Usually safe (these are not detection targets), but document and inspect.               |

`overall_pass = True` means the guide cleared all three screens at the
configured thresholds; it is still only an in silico candidate.

---

## 7. Folding results into the preprint (Figure 2)

Suggested figure layout for the framework preprint:

- **Panel A.** Stacked bar of `n_hits_ge18` per guide, coloured by database.
- **Panel B.** Heatmap of best %identity per (guide x database).
- **Panel C.** Seed-engagement matrix (boolean), one row per guide.
- **Panel D.** Summary donut: `overall_pass` count vs. `fail` count, with
  per-database breakdown.

Source data: `results/stress_test_summary.tsv` (panel A, B, D) and
`results/stress_test_<timestamp>.tsv` (panel C, plus supplementary table).

In the Methods section, describe:
1. Which references and accessions were used (Section 2 of this README).
2. BLAST+ version and `blastn-short` parameters
   (`word_size=7`, `evalue=10`, `dust=no`).
3. The three thresholds (Section 6).
4. Restate that these are **in silico** screens and that wet-lab
   cross-reactivity panels remain required.
