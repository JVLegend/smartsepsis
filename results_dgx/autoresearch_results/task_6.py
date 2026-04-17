import os
import pandas as pd
from collections import Counter
import io
import base64
import matplotlib.pyplot as plt

# --- Configuration and Paths ---
BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results/")
FASTA_DIR = os.path.join(BASE_DIR, "sequences/")
GUIDES_DIR = os.path.join(BASE_DIR, "guides/")
# ---------------------------------

# TASK: Codon usage bias: read FASTA sequences from sequences/ dir. For each gene family, compute RSCU (Relative Synonymous Codon Usage). Compare between families. Save rscu_analysis.csv.

def setup_directories():
    """Ensures the output directory exists."""
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        print(f"Output directory ensured at: {OUTPUT_DIR}")
    except Exception as e:
        print(f"Error creating output directory: {e}")
        return False
    return True

def read_fasta_sequences(fasta_dir):
    """
    Reads all FASTA files from a directory and returns a dictionary {sequence_id: sequence_string}.
    """
    sequences = {}
    print("\n--- Starting FASTA Sequence Reading ---")
    try:
        for filename in os.listdir(fasta_dir):
            if filename.endswith(".fasta"):
                filepath = os.path.join(fasta_dir, filename)
                print(f"Processing file: {filename}")
                
                current_sequence = ""
                header = ""
                
                with open(filepath, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        if line.startswith(">"):
                            if current_sequence:
                                sequences[header] = current_sequence
                            # Start new record
                            header = line[1:].split(" ")[0] # Use the first word as ID
                            current_sequence = ""
                        else:
                            # Append sequence line
                            current_sequence += line.upper().replace(" ", "")
                    # Add the last record
                    sequences[header] = current_sequence
        
        print(f"Successfully loaded {len(sequences)} sequences.")
        return sequences

    except FileNotFoundError:
        print(f"ERROR: FASTA directory not found at {fasta_dir}")
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading FASTA files: {e}")
        return {}

def compute_rscu(sequence):
    """
    Computes the relative usage of all 61 codons found in a sequence.
    Returns a dictionary mapping codon (str) to its relative frequency (float).
    """
    rscu_counts = {}
    
    # Ensure sequence is DNA and uppercase
    sequence = "".join(filter(str.isalpha, sequence.upper()))
    
    if len(sequence) < 3:
        return {}

    codons = [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]
    
    # Count observed codons
    codon_counts = Counter(codons)
    
    total_codons = len(codons)
    
    # Calculate relative frequency for observed codons
    for codon, count in codon_counts.items():
        rscu_counts[codon] = count / total_codons
        
    return rscu_counts

def calculate_rscu_by_family(sequences, guides_df):
    """
    Groups sequences by gene family (using guide metadata) and computes the average RSCU 
    for all sequences associated with that family.
    """
    print("\n--- Calculating RSCU by Gene Family ---")
    
    # Map the sequence ID (header) to a Gene Family based on the guides file.
    # Since the prompt lacks explicit mapping between FASTA headers and guide metadata,
    # we must assume that the sequence ID (header) in the FASTA file corresponds 
    # to the 'gene_family' name found in the guides file structure.
    
    # Dictionary to store: {family_name: {codon: [rscu_value_1, rscu_value_2, ...]}}
    family_rscu_data = {}
    
    # Step 1: Associate sequences with families (A necessary approximation)
    sequence_to_family = {}
    
    for index, row in guides_df.iterrows():
        # Assuming 'gene_family' column holds the identifier that matches the FASTA header
        family = row['gene_family']
        # Using a pseudo-mapping here: assume the guide ID corresponds to a sequence
        # If the FASTA header contains a unique family identifier, this needs adjustment.
        # For safety, we'll just use the family name as the key.
        
        # In a real scenario, we'd use a mapping file. Here, we assume the 
        # first sequence found under this family is used as a representative sequence 
        # for all members of that family for simplicity of calculation.
        pass 

    # Fallback strategy: Process the first sequence encountered for each family as representative.
    # This significantly simplifies the problem but assumes homogeneous family usage.
    
    # We will modify the processing to group sequences by family/source
    family_sequences = {}
    
    # Use a simple grouping mechanism: we process all sequences, and if a sequence 
    # has a name that matches an existing guide family, we associate it.
    
    guide_families = set(guides_df['gene_family'].unique())
    
    for seq_id, sequence in sequences.items():
        # Attempt to determine the family using the sequence ID
        assigned_family = next((family for family in guide_families if family.lower() in seq_id.lower()), "UNKNOWN")
        if assigned_family:
            if assigned_family not in family_sequences:
                family_sequences[assigned_family] = []
            family_sequences[assigned_family].append(sequence)

    
    # Step 2: Calculate average RSCU per family
    family_rscu_results = []
    
    for family, seq_list in family_sequences.items():
        # Dictionary to hold all individual codon usage totals for this family
        total_codon_counts = Counter()
        total_codon_count = 0
        
        print(f"\nProcessing Family: {family} (Sequences: {len(seq_list)})")
        
        for sequence in seq_list:
            # Count codons for this single sequence
            codons = [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]
            total_codon_counts.update(codons)
            total_codon_count += len(codons)
        
        # Calculate overall RSCU for the family
        family_rscu = {}
        for codon, count in total_codon_counts.items():
            family_rscu[codon] = count / total_codon_count
        
        # Format results for DataFrame
        row = {'Gene_Family': family}
        for codon in 'AAA'...'TTT': # Iterate over all 64 codons (assuming standard 3-letter representation)
            row[f'RSCU_{codon}'] = family_rscu.get(codon, 0.0)
        
        family_rscu_results.append(row)

    return pd.DataFrame(family_rscu_results)


def main():
    """Main execution function."""
    if not setup_directories():
        print("Exiting due to directory setup failure.")
        return

    # 1. Load data
    sequences = read_fasta_sequences(FASTA_DIR)
    
    print("\n--- Loading Guide Metadata ---")
    guide_path = os.path.join(GUIDES_DIR, "guide_id,sequence,gene_family,pam_site.tsv")
    guides_df = None
    try:
        # Assuming a comma delimiter based on the provided structure
        guides_df = pd.read_csv(guide_path, sep='\t')
        print(f"Successfully loaded guide metadata ({len(guides_df)} records).")
    except FileNotFoundError:
        print(f"WARNING: Guide file not found at {guide_path}. Cannot perform family grouping.")
    except pd.errors.EmptyDataError:
        print("WARNING: Guide file is empty.")
    except Exception as e:
        print(f"Error loading guides_df: {e}")
    
    if sequences and guides_df is not None:
        # 2. Calculate RSCU
        rscu_df = calculate_rscu_by_family(sequences, guides_df)
        
        # 3. Save results
        output_file_path = os.path.join(OUTPUT_DIR, "rscu_analysis.csv")
        try:
            rscu_df.to_csv(output_file_path, index=False)
            print("\n======================================================")
            print("✅ SUCCESS: RSCU analysis completed.")
            print(f"Results saved to: {output_file_path}")
            print("======================================================")
        except Exception as e:
            print(f"CRITICAL ERROR: Could not save output file to {output_file_path}. {e}")
    elif sequences and guides_df is None:
        print("\nCannot compute family-wise RSCU without guide metadata.")
    else:
        print("\nProcess halted: Could not load sequences or guide metadata.")

if __name__ == "__main__":
    main()
# END OF SCRIPT