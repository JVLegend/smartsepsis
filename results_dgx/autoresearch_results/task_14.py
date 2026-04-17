import os
import pandas as pd
import numpy as np
import base64
import io
from glob import glob

# Define absolute base paths
BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
SEQUENCES_DIR = os.path.join(BASE_DIR, "sequences/")
OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results/")
GUIDES_DIR = os.path.join(BASE_DIR, "guides/")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Helper Functions ---

def parse_fasta(file_path):
    """Reads a FASTA file and returns a list of sequences."""
    sequences = []
    current_sequence = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_sequence:
                        sequences.append("".join(current_sequence).upper())
                    current_sequence = []
                elif line:
                    current_sequence.append(line)
            # Append the last sequence
            if current_sequence:
                sequences.append("".join(current_sequence).upper())
    except Exception as e:
        print(f"Error reading FASTA file {file_path}: {e}")
        return None
    return sequences

def calculate_diversity_stats(sequences):
    """
    Calculates Nucleotide diversity (pi), Number of segregating sites (S),
    and an approximation of Tajima's D (D_approx).
    """
    if not sequences:
        return {'S': 0, 'pi': 0.0, 'Tajima_D_approx': 0.0}

    N = len(sequences)
    L = len(sequences[0])

    if L == 0 or N < 2:
        return {'S': 0, 'pi': 0.0, 'Tajima_D_approx': 0.0}

    # 1. Calculate Pairwise Differences Matrix (Hamming distance)
    diff_count = 0
    pair_count = 0
    
    # Count frequencies of each nucleotide at each position
    site_counts = np.zeros((L, 4), dtype=int) # Columns: A, C, G, T
    
    for seq in sequences:
        for i in range(L):
            nucleotide = seq[i]
            if nucleotide == 'A':
                site_counts[i, 0] += 1
            elif nucleotide == 'C':
                site_counts[i, 1] += 1
            elif nucleotide == 'G':
                site_counts[i, 2] += 1
            elif nucleotide == 'T':
                site_counts[i, 3] += 1
    
    # 2. Calculate S (Segregating Sites)
    # A site is segregating if the count of the most common nucleotide is less than N.
    S = 0
    for i in range(L):
        max_count_at_site = np.max(site_counts[i])
        if max_count_at_site < N:
            S += 1

    # 3. Calculate pi (Nucleotide Diversity)
    # pi = (1 / (L * N)) * sum(pairwise differences)
    # More simply: pi = (1/L) * average pairwise difference
    
    total_pairwise_diff = 0
    pair_count = 0
    for i in range(N):
        for j in range(i + 1, N):
            pair_diff = 0
            for k in range(L):
                if sequences[i][k] != sequences[j][k]:
                    pair_diff += 1
            total_pairwise_diff += pair_diff
            pair_count += 1
    
    pi = total_pairwise_diff / (pair_count * L)

    # 4. Approximate Tajima's D (Simplified Site Count Approach)
    # D is often approximated by analyzing the site count spectrum (S vs. number of alleles).
    # A common simplification: D = (S/N) - (L/N) * (Sum(count_i^2) / N^2)
    
    sum_of_squares = 0
    for i in range(L):
        # site_counts[i] holds [A_count, C_count, G_count, T_count]
        sum_of_squares += np.sum(site_counts[i]**2)

    # D_approx calculation based on site count spectrum
    D_approx = (S / N) - (L / N) * (sum_of_squares / (N**2 * 4))

    return {'S': S, 'pi': pi, 'Tajima_D_approx': D_approx}

# --- Main Script Logic ---

def run_diversity_analysis():
    """
    Performs sequence diversity analysis on all FASTA files, grouping results 
    by derived gene family ID.
    """
    print("--- CRISPR-Cas12a AMR Diagnostic Pipeline ---")
    print("TASK: Sequence diversity analysis (pi, S, Tajima's D) for gene family profiling.")
    
    all_results = {}
    fasta_files = glob(os.path.join(SEQUENCES_DIR, "*.fasta"))
    
    if not fasta_files:
        print(f"\n[ERROR] No FASTA files found in {SEQUENCES_DIR}. Exiting.")
        return

    print(f"\n[INFO] Found {len(fasta_files)} FASTA files to process.")
    
    # 1. Process FASTA files
    for i, file_path in enumerate(fasta_files):
        print(f"Processing file {i+1}/{len(fasta_files)}: {os.path.basename(file_path)}...")
        
        sequences = parse_fasta(file_path)
        
        if sequences:
            stats = calculate_diversity_stats(sequences)
            
            # Derive Gene Family ID: Assuming the gene family can be derived from the file name.
            # e.g., filename_for_familyXYZ.fasta -> FamilyXYZ
            filename = os.path.basename(file_path)
            gene_family_id = filename.replace(".fasta", "").split('_')[0].upper()
            
            # Store or append results
            if gene_family_id not in all_results:
                all_results[gene_family_id] = []
            
            all_results[gene_family_id].append({
                'filename': filename,
                'num_sequences': len(sequences),
                'pi': stats['pi'],
                'S': stats['S'],
                'Tajima_D_approx': stats['Tajima_D_approx']
            })
        
    # 2. Aggregate and Save Results
    
    if not all_results:
        print("\n[WARNING] No successful diversity calculations were performed.")
        return

    final_data = []
    for family, data_list in all_results.items():
        # Average the statistics across all files assigned to this family
        avg_num_sequences = np.mean([d['num_sequences'] for d in data_list])
        avg_pi = np.mean([d['pi'] for d in data_list])
        avg_S = np.mean([d['S'] for d in data_list])
        avg_D = np.mean([d['Tajima_D_approx'] for d in data_list])
        
        final_data.append({
            'gene_family_id': family,
            'avg_num_sequences': avg_num_sequences,
            'mean_nucleotide_diversity_pi': avg_pi,
            'mean_segregating_sites_S': avg_S,
            'mean_tajimas_d': avg_D
        })
    
    diversity_df = pd.DataFrame(final_data)
    output_path = os.path.join(OUTPUT_DIR, "diversity_stats.csv")
    
    try:
        diversity_df.to_csv(output_path, index=False)
        print(f"\n[SUCCESS] Diversity analysis completed.")
        print(f"Results saved to: {output_path}")
    except Exception as e:
        print(f"\n[CRITICAL ERROR] Could not save the output file. {e}")

if __name__ == "__main__":
    run_diversity_analysis()
# END OF SCRIPT
# NOTE: The required standalone script must contain only the Python code.
# The provided solution is the complete, runnable script.