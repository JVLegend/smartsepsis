import pandas as pd
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

# Define the base directory
BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
GUIDES_DIR = os.path.join(BASE_DIR, "guides/")
SEQUENCES_DIR = os.path.join(BASE_DIR, "sequences/")
OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results/")

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def read_all_fasta_sequences(directory):
    """
    Reads all FASTA files in the given directory and maps sequence to its conceptual gene/variant.
    Since the FASTA files represent individual gene variants, we'll extract metadata 
    (e.g., filename or assuming the file name dictates the 'gene' or 'variant').
    """
    print(f"--- Reading FASTA sequences from {directory} ---")
    target_sequences = []
    
    for filepath in glob.glob(os.path.join(directory, "*.fasta")):
        print(f"Processing FASTA file: {os.path.basename(filepath)}")
        try:
            # Assuming each FASTA file contains at least one sequence
            records = list(SeqIO.parse(filepath, "fasta"))
            if records:
                seq_record = records[0]
                target_sequences.append({
                    'variant': os.path.splitext(os.path.basename(filepath))[0],
                    'sequence': str(seq_record.seq)
                })
        except Exception as e:
            print(f"Error reading FASTA file {filepath}: {e}")
            continue
    
    return pd.DataFrame(target_sequences)

def calculate_mismatches(seq1, seq2, max_mismatches=3):
    """
    Calculates the Hamming distance (mismatches) between two sequences of equal length.
    If lengths differ, we pad or truncate to the minimum length for comparison.
    """
    len1 = len(seq1)
    len2 = len(seq2)
    
    min_len = min(len1, len2)
    
    if min_len == 0:
        return 0, 0 # No meaningful comparison
        
    # Compare up to the minimum length
    mismatches = 0
    for i in range(min_len):
        if seq1[i] != seq2[i]:
            mismatches += 1
            
    # Check if the mismatch count exceeds the limit
    if mismatches > max_mismatches:
        return mismatches, False
    
    return mismatches, True


def analyze_cross_reactivity(guide_df, target_df, max_mismatches=3):
    """
    Analyzes off-target hits for each guide against all target sequences, 
    excluding hits originating from the same gene family as the guide.
    """
    cross_reactivity_hits = []
    
    print("\n--- Starting Cross-Reactivity Analysis ---")
    
    # Since the FASTA files don't explicitly contain gene family info, 
    # we must assume that the target sequences are variants derived from 
    # distinct genes, and we need a proxy for "other gene families". 
    # For this exercise, we will assume the 'variant' name (derived from the filename) 
    # acts as a proxy for the gene family/source identifier.
    
    # Group 1: Guides (source: guide_df)
    # Group 2: Targets (source: target_df)

    for index, guide_row in guide_df.iterrows():
        spacer_seq = guide_row['spacer_seq']
        guide_gene = guide_row['gene']
        
        print(f"\n[Processing Guide from Gene: {guide_gene}, Spacer: {spacer_seq}]")
        
        for index, target_row in target_df.iterrows():
            target_variant = target_row['variant']
            target_seq = target_row['sequence']
            
            # CRITICAL RULE CHECK: Must compare against OTHER gene families.
            # We assume the target gene source is encapsulated by the 'variant' name.
            # For simplicity, we require the source (variant name) to not match 
            # a known gene source (guide_gene). This is a heuristic based on available data.
            # A more robust solution would require explicit gene family linkage.
            if target_variant.startswith(guide_gene) or target_variant == guide_gene:
                # This heuristic check prevents self-cross-reactivity within the same apparent gene family
                continue 
            
            # Calculate mismatches
            mismatches, is_hit = calculate_mismatches(
                spacer_seq, 
                target_seq, 
                max_mismatches=max_mismatches
            )
            
            if is_hit:
                hit_info = {
                    'guide_gene': guide_row['gene'],
                    'guide_spacer_seq': spacer_seq,
                    'target_variant': target_variant,
                    'target_sequence': target_seq,
                    'mismatches': mismatches,
                    'is_cross_reactive': True
                }
                cross_reactivity_hits.append(hit_info)

    return pd.DataFrame(cross_reactivity_hits)

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    
    # 1. Load Guide Data
    print("Loading guide sequence data...")
    guide_files = glob.glob(os.path.join(GUIDES_DIR, "*.tsv"))
    guide_list = []
    
    for filepath in guide_files:
        try:
            # Use sep='\t' and specific required columns
            df = pd.read_csv(filepath, sep='\t')
            print(f"Successfully loaded guide file: {os.path.basename(filepath)}")
            print(f'Columns: {df.columns.tolist()}')
            
            # Assuming we concatenate all guide data into one master guide DataFrame
            guide_list.append(df[['gene', 'spacer_seq']].copy())
            
        except FileNotFoundError:
            print(f"ERROR: File not found at {filepath}")
        except Exception as e:
            print(f"ERROR processing {filepath}: {e}")
            
    if not guide_list:
        print("FATAL ERROR: No guide data was successfully loaded. Exiting.")
        exit()
        
    master_guide_df = pd.concat(guide_list, ignore_index=True)
    print(f"\nTotal guide entries loaded: {len(master_guide_df)}")

    # 2. Load Target Sequences (FASTA files)
    target_df = read_all_fasta_sequences(SEQUENCES_DIR)
    
    if target_df.empty:
        print("\nWARNING: No target sequences were loaded from FASTA files. Cannot perform cross-reactivity analysis.")
        exit()

    # 3. Run Analysis
    try:
        cross_reactivity_df = analyze_cross_reactivity(
            master_guide_df, 
            target_df, 
            max_mismatches=3
        )
    except Exception as e:
        print(f"\nCRITICAL ERROR during cross-reactivity analysis: {e}")
        cross_reactivity_df = pd.DataFrame()

    # 4. Save Results
    output_filename = os.path.join(OUTPUT_DIR, "cross_reactivity.csv")
    
    if not cross_reactivity_df.empty:
        cross_reactivity_df.to_csv(output_filename, index=False)
        print("\n" + "="*60)
        print(f"✅ Success! Cross-reactivity analysis complete.")
        print(f"Results saved to: {output_filename}")
        print(f"Total hits found: {len(cross_reactivity_df)}")
        print("="*60)
    else:
        print("\n" + "="*60)
        print("⚠️ Warning: No cross-reactive hits found or an error occurred.")
        print("No output file was generated.")
        print("="*60)
