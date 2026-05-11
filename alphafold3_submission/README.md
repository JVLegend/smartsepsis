# AlphaFold 3 Server — Ternary Complex Submission Pack

**Goal:** model **LbCas12a + crRNA + dsDNA target** ternary complex for the top
mecA guide candidates, to assess structural feasibility before wet-lab work.
This becomes Figure 3 of `preprint/draft_v0.md`.

## Status

| File | What it is |
|---|---|
| `cas12a_lba_protein.fasta` | LbCas12a (UniProt A0Q7Q2), 1228 aa — paste in AF3 "Protein" chain |
| `top5_mecA_crRNAs.fasta` | 5 crRNAs (scaffold + 20-nt spacer) — paste in "RNA" chain |
| `targets/target_dsDNA_*.fasta` | Top + bottom strands of target dsDNA — paste in "DNA" chain |
| `targets/SUMMARY.tsv` | PAM regions and positions confirmed against NG_047945.1 (mecA RefSeq) |

## Submission ready — 3 ternary jobs

For these 3 guides, target dsDNA was extracted from the mecA RefSeq
(NG_047945.1) at the exact positions identified by `design_guides.py`:

| guide_id | PAM (TTTV) | strand | position |
|---|---|---|---|
| **mecA-0** | ATTTTA → TTTA | + | 465 |
| **mecA-1** | CATTTA → TTTA | + | 1114 |
| **mecA-2** | TATTTC → TTTC | − | 541 |

All three have canonical Cas12a PAM `TTTV` (V = A/C/G) immediately 5' of the
protospacer on the non-spacer strand. **Submit these first.**

## Pending — mecA1-0 and mecA2-0

mecA1 and mecA2 are paralogs/homologs to canonical mecA with distinct
sequences. Their target regions are in different RefSeq entries (full
*Staphylococcus* contigs, not isolated gene references). Defer for a second
batch — fetch the exact spacer position in each respective genome,
extract ±15 nt flanks, and submit separately.

## How to submit (manual, AlphaFold Server)

1. Open https://alphafoldserver.com (Google account required).
2. Click **"New job"** → set **3 chains**:
   - **Protein**: paste content of `cas12a_lba_protein.fasta` (one chain, copies=1)
   - **RNA**: paste the relevant `>crRNA_mecA-X` sequence from `top5_mecA_crRNAs.fasta`
   - **DNA**: paste BOTH lines (`target_top_...` and `target_bot_...`) from `targets/target_dsDNA_mecA-X.fasta` — 2 chains of DNA
3. Use default seed; submit. Each job ~10 min.
4. Download CIF + JSON; save to `alphafold3_submission/results/<guide_id>/`.

## After AF3 returns

For each complex, manually inspect or script-extract:
- pLDDT mean and ipTM (overall + per-chain)
- crRNA–target dsDNA base-pair coverage (visual check of duplex formation)
- Cas12a active-site geometry near the scissile phosphate of the non-target strand

Use ChimeraX or PyMOL for visualization. Top-pLDDT complex becomes Figure 3
of the preprint with a 4-panel layout: structure cartoon, pLDDT heatmap,
RNA-DNA duplex zoom, active-site geometry.

## In silico only

These predictions do not confirm cleavage. They confirm structural
plausibility of the designed complex. Experimental validation remains
pending — see `CONTRIBUTING.md` for partnership inquiries.
