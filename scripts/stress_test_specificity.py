#!/usr/bin/env python3
"""
stress_test_specificity.py
==========================

In silico off-target / specificity stress-test for CRISPR-Cas12a guide RNAs
designed by the SmartSepsis-Oph pipeline.

The script queries each guide spacer against three locally-built BLAST+
nucleotide databases:

  (a) human_genome              - GRCh38, to flag host DNA off-targets.
  (b) ocular_microbiome         - 16S / WGS of commensal eye bacteria.
  (c) ocular_fungi_acanthamoeba - Acanthamoeba + ocular fungal references
                                   (Aspergillus fumigatus, Candida albicans,
                                   Fusarium solani) used as negative controls.

For every (guide, database) pair we compute:
  - number of hits with >= 18 nt of contiguous complementarity
  - best-hit identity %
  - whether the seed region (spacer positions 1-8 for Cas12a) is involved
    in the alignment (a hit that engages the seed is a higher concern than
    a hit that only matches the distal end of the spacer).

Two TSV reports are produced under results/:

  results/stress_test_<timestamp>.tsv
      one row per (guide, database, hit). Empty hits are recorded too.

  results/stress_test_summary.tsv
      one row per guide, aggregating pass / fail per database under the
      thresholds defined in DEFAULT_THRESHOLDS (see below).

IMPORTANT
---------
All results are **in silico**. This script does NOT establish clinical or
analytical specificity on its own. A negative BLAST result is a necessary
but not sufficient condition; wet-lab cross-reactivity panels are still
required before any clinical claim.

Usage
-----
    python scripts/stress_test_specificity.py \\
        --input data/guides_panel.csv \\
        --db-human   /path/to/blast_dbs/human_genome \\
        --db-microb  /path/to/blast_dbs/ocular_microbiome \\
        --db-fungi   /path/to/blast_dbs/ocular_fungi_acanthamoeba \\
        --outdir results/

Input file (CSV or parquet) must contain the columns:
    gene_family, variant, spacer_sequence, strand
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

# --------------------------------------------------------------------------
# Thresholds for pass/fail. Conservative defaults; tune per use case.
# --------------------------------------------------------------------------
DEFAULT_THRESHOLDS = {
    # A guide FAILS for a database if ANY of the following is true:
    "max_hits_18nt": 0,        # any hit with >=18 nt complementarity is a fail
    "max_identity_pct": 85.0,  # best-hit identity must be < this %
    "seed_strict": True,       # any seed-region (pos 1-8) hit is an automatic fail
}

# Cas12a "seed" region. Cas12a is most sensitive to mismatches in the
# PAM-proximal 1-8 nt of the spacer. We use this for the seed-engagement flag.
SEED_START_1BASED = 1
SEED_END_1BASED = 8

# Deterministic seed for any stochastic ops (BLAST itself is deterministic
# given the same database + parameters, but we set PYTHONHASHSEED upstream).
os.environ.setdefault("PYTHONHASHSEED", "0")


# --------------------------------------------------------------------------
# Data classes
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Guide:
    gene_family: str
    variant: str
    spacer_sequence: str
    strand: str

    @property
    def guide_id(self) -> str:
        return f"{self.gene_family}|{self.variant}|{self.strand}"


@dataclass
class BlastHit:
    guide_id: str
    database: str
    subject_id: str
    subject_title: str
    pident: float           # % identity over the aligned segment
    length: int             # alignment length (nt)
    mismatch: int
    gapopen: int
    qstart: int             # 1-based on spacer
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    engages_seed: bool


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------
def load_guides(path: Path) -> list[Guide]:
    """Load guides from a CSV or parquet file."""
    if not path.exists():
        raise FileNotFoundError(path)

    suffix = path.suffix.lower()
    rows: list[dict]
    if suffix in (".csv", ".tsv"):
        delim = "\t" if suffix == ".tsv" else ","
        with path.open() as f:
            rows = list(csv.DictReader(f, delimiter=delim))
    elif suffix in (".parquet", ".pq"):
        try:
            import pandas as pd  # local import; only needed for parquet
        except ImportError as exc:
            raise RuntimeError(
                "pandas is required to read parquet input"
            ) from exc
        rows = pd.read_parquet(path).to_dict(orient="records")
    else:
        raise ValueError(f"Unsupported input extension: {suffix}")

    required = {"gene_family", "variant", "spacer_sequence", "strand"}
    if rows and not required.issubset(rows[0].keys()):
        missing = required - set(rows[0].keys())
        raise ValueError(f"Input missing required columns: {missing}")

    guides: list[Guide] = []
    for r in rows:
        spacer = str(r["spacer_sequence"]).strip().upper()
        if not spacer or any(b not in "ACGTN" for b in spacer):
            logging.warning(
                "Skipping malformed spacer for %s|%s: %r",
                r.get("gene_family"), r.get("variant"), spacer,
            )
            continue
        if not (20 <= len(spacer) <= 25):
            logging.warning(
                "Spacer length %d outside [20,25] for %s|%s",
                len(spacer), r.get("gene_family"), r.get("variant"),
            )
        guides.append(Guide(
            gene_family=str(r["gene_family"]),
            variant=str(r["variant"]),
            spacer_sequence=spacer,
            strand=str(r["strand"]),
        ))
    return guides


def write_query_fasta(guides: Iterable[Guide], out: Path) -> None:
    """Write all spacers to a single FASTA for one BLAST run per DB."""
    with out.open("w") as f:
        for g in guides:
            f.write(f">{g.guide_id}\n{g.spacer_sequence}\n")


# --------------------------------------------------------------------------
# BLAST
# --------------------------------------------------------------------------
def run_blastn(query_fa: Path, db: str, out_tsv: Path,
               word_size: int = 7, evalue: float = 10.0,
               threads: int = 4) -> None:
    """
    Invoke blastn (BLAST+ command line) against a local nucleotide DB.

    word_size=7 is short enough to surface partial matches relevant to a
    20-25 nt spacer; evalue is loose because we want to capture marginal
    off-targets too, then filter in Python.
    """
    if shutil.which("blastn") is None:
        raise RuntimeError(
            "blastn (NCBI BLAST+) not found in PATH. "
            "Install BLAST+ and rebuild databases per README_stress_test.md."
        )

    cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", str(query_fa),
        "-db", db,
        "-out", str(out_tsv),
        "-outfmt",
        "6 qseqid sseqid stitle pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore",
        "-word_size", str(word_size),
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-dust", "no",
    ]
    logging.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_blast_tsv(path: Path, database_label: str) -> list[BlastHit]:
    """Parse blastn outfmt-6 with stitle (13 cols)."""
    hits: list[BlastHit] = []
    if not path.exists() or path.stat().st_size == 0:
        return hits
    with path.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 13:
                continue
            (qseqid, sseqid, stitle, pident, length, mismatch, gapopen,
             qstart, qend, sstart, send, evalue, bitscore) = parts[:13]
            qs, qe = int(qstart), int(qend)
            engages_seed = not (qe < SEED_START_1BASED or qs > SEED_END_1BASED)
            hits.append(BlastHit(
                guide_id=qseqid,
                database=database_label,
                subject_id=sseqid,
                subject_title=stitle,
                pident=float(pident),
                length=int(length),
                mismatch=int(mismatch),
                gapopen=int(gapopen),
                qstart=qs,
                qend=qe,
                sstart=int(sstart),
                send=int(send),
                evalue=float(evalue),
                bitscore=float(bitscore),
                engages_seed=engages_seed,
            ))
    return hits


# --------------------------------------------------------------------------
# Aggregation / pass-fail
# --------------------------------------------------------------------------
def summarise(guides: list[Guide],
              hits_by_db: dict[str, list[BlastHit]],
              thresholds: dict) -> list[dict]:
    """Per-guide x per-database pass/fail summary."""
    rows: list[dict] = []
    for g in guides:
        row: dict = {
            "guide_id": g.guide_id,
            "gene_family": g.gene_family,
            "variant": g.variant,
            "strand": g.strand,
            "spacer_sequence": g.spacer_sequence,
        }
        any_fail = False
        for db_label, hits in hits_by_db.items():
            ghits = [h for h in hits if h.guide_id == g.guide_id]
            ghits_18 = [h for h in ghits if h.length >= 18]
            best_pident = max((h.pident for h in ghits), default=0.0)
            seed_hits = [h for h in ghits_18 if h.engages_seed]

            fail = (
                len(ghits_18) > thresholds["max_hits_18nt"]
                or best_pident >= thresholds["max_identity_pct"]
                or (thresholds["seed_strict"] and len(seed_hits) > 0)
            )
            any_fail = any_fail or fail

            row[f"{db_label}_n_hits_ge18"] = len(ghits_18)
            row[f"{db_label}_best_pident"] = round(best_pident, 2)
            row[f"{db_label}_seed_hits"] = len(seed_hits)
            row[f"{db_label}_pass"] = (not fail)

        row["overall_pass"] = (not any_fail)
        rows.append(row)
    return rows


def write_hits_tsv(path: Path, guides: list[Guide],
                   hits_by_db: dict[str, list[BlastHit]]) -> None:
    """One row per (guide, database, hit). Empty hits get a sentinel row."""
    header = [
        "guide_id", "gene_family", "variant", "strand", "spacer_sequence",
        "database", "subject_id", "subject_title", "pident", "length",
        "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore", "engages_seed",
    ]
    guide_lookup = {g.guide_id: g for g in guides}
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for db_label, hits in hits_by_db.items():
            seen = set()
            for h in hits:
                g = guide_lookup.get(h.guide_id)
                if g is None:
                    continue
                seen.add(h.guide_id)
                w.writerow([
                    g.guide_id, g.gene_family, g.variant, g.strand,
                    g.spacer_sequence, db_label, h.subject_id, h.subject_title,
                    f"{h.pident:.2f}", h.length, h.mismatch, h.gapopen,
                    h.qstart, h.qend, h.sstart, h.send,
                    f"{h.evalue:.2e}", f"{h.bitscore:.2f}", h.engages_seed,
                ])
            # Sentinel "no hit" rows so the report is complete.
            for g in guides:
                if g.guide_id not in seen:
                    w.writerow([
                        g.guide_id, g.gene_family, g.variant, g.strand,
                        g.spacer_sequence, db_label, "", "NO_HIT",
                        "0.00", 0, 0, 0, 0, 0, 0, 0, "", "0.00", False,
                    ])


def write_summary_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    header = list(rows[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="In silico off-target stress-test for Cas12a guides.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", required=True, type=Path,
                   help="CSV or parquet of guides "
                        "(cols: gene_family, variant, spacer_sequence, strand).")
    p.add_argument("--db-human", required=True,
                   help="Path/prefix of human_genome BLAST DB.")
    p.add_argument("--db-microb", required=True,
                   help="Path/prefix of ocular_microbiome BLAST DB.")
    p.add_argument("--db-fungi", required=True,
                   help="Path/prefix of ocular_fungi_acanthamoeba BLAST DB.")
    p.add_argument("--outdir", default=Path("results"), type=Path)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--word-size", type=int, default=7)
    p.add_argument("--evalue", type=float, default=10.0)
    p.add_argument("--max-hits-18nt", type=int,
                   default=DEFAULT_THRESHOLDS["max_hits_18nt"])
    p.add_argument("--max-identity-pct", type=float,
                   default=DEFAULT_THRESHOLDS["max_identity_pct"])
    p.add_argument("--seed-strict", action="store_true", default=True)
    p.add_argument("--no-seed-strict", dest="seed_strict",
                   action="store_false")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    logging.info("Loading guides from %s", args.input)
    guides = load_guides(args.input)
    logging.info("Loaded %d guides", len(guides))
    if not guides:
        logging.error("No guides to test. Aborting.")
        return 2

    # Write one combined FASTA and run blastn once per DB.
    work = args.outdir / f".tmp_{timestamp}"
    work.mkdir(parents=True, exist_ok=True)
    query_fa = work / "guides.fasta"
    write_query_fasta(guides, query_fa)

    db_map = {
        "human_genome": args.db_human,
        "ocular_microbiome": args.db_microb,
        "ocular_fungi_acanthamoeba": args.db_fungi,
    }

    hits_by_db: dict[str, list[BlastHit]] = {}
    for label, db_path in db_map.items():
        out_tsv = work / f"blast_{label}.tsv"
        try:
            run_blastn(query_fa, db_path, out_tsv,
                       word_size=args.word_size, evalue=args.evalue,
                       threads=args.threads)
        except subprocess.CalledProcessError as exc:
            logging.error("blastn failed for %s: %s", label, exc)
            return 3
        hits_by_db[label] = parse_blast_tsv(out_tsv, label)
        logging.info("%s: %d raw hits", label, len(hits_by_db[label]))

    thresholds = {
        "max_hits_18nt": args.max_hits_18nt,
        "max_identity_pct": args.max_identity_pct,
        "seed_strict": args.seed_strict,
    }

    hits_path = args.outdir / f"stress_test_{timestamp}.tsv"
    summary_path = args.outdir / "stress_test_summary.tsv"
    write_hits_tsv(hits_path, guides, hits_by_db)
    summary_rows = summarise(guides, hits_by_db, thresholds)
    write_summary_tsv(summary_path, summary_rows)

    n_pass = sum(1 for r in summary_rows if r["overall_pass"])
    logging.info("Done. %d / %d guides PASS all three databases.",
                 n_pass, len(summary_rows))
    logging.info("Hit-level report:  %s", hits_path)
    logging.info("Per-guide summary: %s", summary_path)
    logging.info("NOTE: results are in silico only; wet-lab cross-reactivity "
                 "panels are still required.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
