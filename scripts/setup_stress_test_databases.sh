#!/usr/bin/env bash
# SmartSepsis-Oph — stress-test database setup
#
# Downloads reference databases and builds local BLAST+ indices for the
# specificity stress-test (scripts/stress_test_specificity.py).
#
# Total disk required: ~5-8 GB (GRCh38 dominates).
# Time: 30 min to several hours depending on bandwidth.
#
# Requirements:
#   - BLAST+ (makeblastdb, blastn)  — brew install blast OR apt install ncbi-blast+
#   - datasets CLI (NCBI)           — https://github.com/ncbi/datasets
#   - curl, gunzip
#
# Usage:
#   bash scripts/setup_stress_test_databases.sh /path/to/blast_dbs

set -euo pipefail
BLAST_DIR="${1:-./blast_dbs}"
mkdir -p "$BLAST_DIR" && cd "$BLAST_DIR"
echo "Setting up BLAST databases in: $(pwd)"

# 1. Human genome GRCh38
if [ ! -f human_genome.nhr ]; then
  echo "[1/3] Downloading GRCh38..."
  curl -L -o human.fna.gz \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
  gunzip -k human.fna.gz
  makeblastdb -in human.fna -dbtype nucl -out human_genome -title "GRCh38_p14"
  echo "  human_genome built"
fi

# 2. Ocular microbiome — concatenate commensal eye bacterial refs
if [ ! -f ocular_microbiome.nhr ]; then
  echo "[2/3] Downloading ocular commensal refs (Staphylococcus epidermidis, Cutibacterium acnes, Corynebacterium spp.)..."
  : > ocular_micro.fna
  for ACC in NC_004461.1 NC_006085.1 NZ_CP013128.1; do
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text" >> ocular_micro.fna
    sleep 1
  done
  makeblastdb -in ocular_micro.fna -dbtype nucl -out ocular_microbiome -title "ocular_commensals"
  echo "  ocular_microbiome built. NOTE: This is a minimal set; for full coverage add SILVA 16S + BioProject PRJNA646695 (Doan 2016) reads."
fi

# 3. Acanthamoeba + ocular fungi
if [ ! -f ocular_fungi_acanthamoeba.nhr ]; then
  echo "[3/3] Downloading Acanthamoeba + Aspergillus + Candida + Fusarium..."
  : > fungi_acanth.fna
  for ACC in NW_004457519.1 NC_007194.1 NC_032089.1 NC_037070.1; do
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=text" >> fungi_acanth.fna
    sleep 1
  done
  makeblastdb -in fungi_acanth.fna -dbtype nucl -out ocular_fungi_acanthamoeba -title "ocular_eukaryotes"
  echo "  ocular_fungi_acanthamoeba built"
fi

echo ""
echo "All three databases ready in $(pwd):"
ls -1 *.nhr | sed 's/.nhr//'
echo ""
echo "Next:"
echo "  python scripts/stress_test_specificity.py \\"
echo "    --guides fase7_dgx_results/fase7_results/rnafold_guides.csv \\"
echo "    --blast-dir $(pwd) \\"
echo "    --output-dir results/"
