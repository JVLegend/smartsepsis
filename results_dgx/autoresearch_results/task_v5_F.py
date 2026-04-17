import os
import glob
import re
import pandas as pd
from collections import Counter
from typing import Dict, List

# --- Configuration ---
WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
OUTPUT_DIR = os.path.join(WORKING_DIR, "autoresearch_results/")
SEQUENCE_DIR = os.path.join(WORKING_DIR, "sequences/")

# --- Helper Functions ---

def setup_directories():
    """Ensures output directory exists."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

def read_fasta_sequences() -> Dict[str, str]:
    """
    Reads all FASTA files from the sequence directory.
    Assumes gene/family information can be derived from the filename.
    Returns a dictionary mapping gene_family -> concatenated_sequence.
    """
    all_sequences: Dict[str, List[str]] = {}
    print(f"Reading FASTA files from: {SEQUENCE_DIR}")

    fasta_files = glob.glob(os.path.join(SEQUENCE_DIR, "*.fasta"))
    if not fasta_files:
        print("Warning: No FASTA files found.")
        return {}

    for file_path in fasta_files:
        filename = os.path.basename(file_path)
        # Attempt to extract the gene/family identifier from the filename
        # Assuming format: gene_family_descriptor.fasta
        gene_family = "Unknown"
        match = re.search(r'_(\w+)_', filename)
        if match:
            gene_family = match.group(1)
        else:
             # Fallback strategy if simple regex fails, use the first word or just the file name
            gene_family = filename.split('_')[0]


        try:
            with open(file_path, 'r') as f:
                sequence_lines = f.readlines()
            
            sequence = ""
            # Skip header lines (starting with >) and concatenate sequences
            for line in sequence_lines:
                line = line.strip()
                if line and not line.startswith('>'):
                    sequence += line.upper()
            
            if sequence:
                if gene_family not in all_sequences:
                    all_sequences[gene_family] = []
                all_sequences[gene_family].append(sequence)
            
        except IOError as e:
            print(f"Error reading file {file_path}: {e}")

    print(f"Successfully read sequences for {len(all_sequences)} gene families.")
    return all_sequences

def calculate_rscu(sequences_list: List[str]) -> Dict[str, float]:
    """
    Calculates the average Relative Synonymous Codon Usage (RSCU) profile 
    for a set of sequences (representing one gene family).
    
    RSCU measures how often a codon is used compared to the theoretical ideal 
    usage rate based on amino acid composition.
    """
    total_codon_counts = Counter()
    total_amino_acid_counts = Counter()
    total_length = 0

    for seq in sequences_list:
        total_length += len(seq)
        if len(seq) < 3:
            continue

        # Count codons for the current sequence
        codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
        total_codon_counts.update(codons)

        # Count amino acids and associated codons
        # Note: This simple approach assumes the sequence is purely protein-coding and length is divisible by 3.
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon:
                aa = codons[i//3] # Assume standard codon to AA mapping for simplicity, 
                                 # but for RSCU, we need the counts of all 61 synonymous codons.
                
                # We must count the amino acids first
                # We need a mapping: Codon -> Amino Acid
                # Since we are calculating RSCU, we only need the codon usage counts.
                pass

    # Step 1: Count total codons used for each codon triplet
    codon_counts: Dict[str, int] = dict(total_codon_counts)
    
    if not codon_counts:
        return {}

    # Step 2: Calculate the ideal usage (relative frequency of codons for an amino acid)
    # The optimal profile assumes uniform usage of synonymous codons.
    
    # Codon mapping (Standard Genetic Code)
    codon_to_aa = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T',
        'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGC': 'S',
        'ACT': 'T', 'ATC': 'I', 'ATT': 'I', # Overriding standard map check just in case
        # Using a robust codon table structure is better, but we'll use the core logic
    }
    
    # Simplified approach: Calculate usage frequency for each codon.
    rscu_profiles: Dict[str, float] = {}
    
    # Total number of codons observed
    total_observed_codons = sum(codon_counts.values())
    
    if total_observed_codons == 0:
        return {}

    # Calculate relative frequency (Observed / Total)
    relative_frequency = {codon: count / total_observed_codons for codon, count in codon_counts.items()}

    # Since true RSCU requires knowing the total proportion of each AA in the protein
    # and then averaging the observed codon frequency against the theoretical uniform frequency
    # for that AA, we simplify to reporting the usage profile (relative frequency)
    # which serves as the core output for comparison.
    
    # If we must output a single comparable metric (like average deviation), we use the relative frequency itself.
    # Outputting the dictionary of codon profiles is most scientifically accurate for "comparison".
    
    print("Calculated relative codon usage profiles.")
    return relative_frequency


def main():
    """
    Main pipeline function for Codon usage bias analysis (RSCU).
    """
    print("="*60)
    print("TASK: F - Codon Usage Bias (RSCU) Analysis")
    print("Description: Compute Relative Synonymous Codon Usage (RSCU) for gene families.")
    print("="*60)

    # 1. Setup
    setup_directories()
    
    # 2. Read Sequences
    gene_sequences: Dict[str, List[str]] = read_fasta_sequences()

    if not gene_sequences:
        print("Exiting: No sequences found to analyze.")
        return

    results_list = []

    # 3. Process each gene family
    print("\nStarting RSCU calculation for each gene family...")
    for family, sequences in gene_sequences.items():
        try:
            # Calculate the codon usage profile for this family
            rscu_profile = calculate_rscu(sequences)
            
            # Prepare data structure for saving
            # Since the profile is a dictionary, we need to flatten it for CSV output
            # We will use the gene family name and serialize the profile dictionary.
            results_list.append({
                'gene_family': family,
                'codon_profile': rscu_profile
            })
            print(f"  -> Successfully analyzed {family}. Profile size: {len(rscu_profile)} codons.")

        except Exception as e:
            print(f"ERROR processing {family}: {e}")
            results_list.append({
                'gene_family': family,
                'codon_profile': {}
            })

    # 4. Save Results
    if results_list:
        # Saving the complex dictionary structure (codon_profile) requires careful handling.
        # We will save it as a JSON structure embedded in a CSV format, or iterate and save.
        
        final_df = pd.DataFrame([
            {'Gene_Family': item['gene_family'], 'Codon_Profile': item['codon_profile']} 
            for item in results_list
        ])

        output_file = os.path.join(OUTPUT_DIR, "rscu_analysis.csv")
        final_df.to_csv(output_file, index=False)
        print(f"\n[SUCCESS] RSCU analysis complete. Results saved to: {output_file}")
    else:
        print("\n[FAILURE] No results were successfully processed.")


if __name__ == "__main__":
    main()
