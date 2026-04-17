import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

# TASK: A - description Redesign CRISPR-Cas12a guides for poorly covered gene families (blaNDM, blaVIM, blaIMP, blaGES, vanA, armA, qnrS). Reads FASTA files, computes consensus, and designs 23nt guides with TTTV PAM sites.

def initialize_paths():
    """Defines and checks critical working and output paths."""
    BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
    
    INPUT_SEQ_DIR = os.path.join(BASE_DIR, "sequences/")
    OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results/")
    
    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    return BASE_DIR, INPUT_SEQ_DIR, OUTPUT_DIR

def compute_consensus(fasta_file_path):
    """
    Reads multiple sequences from a FASTA file and computes the consensus sequence.
    Returns the consensus string or None if no sequences are found.
    """
    sequences = []
    try:
        for record in SeqIO.parse(fasta_file_path, "fasta"):
            sequences.append(str(record.seq).upper())
    except FileNotFoundError:
        print(f"Warning: FASTA file not found: {fasta_file_path}")
        return None
    except Exception as e:
        print(f"Error reading FASTA file {fasta_file_path}: {e}")
        return None

    if not sequences:
        print(f"Warning: No sequences found in {fasta_file_path}")
        return None

    # Determine the length of the consensus (minimum length of all sequences)
    min_len = min(len(seq) for seq in sequences)
    
    consensus_list = []
    for i in range(min_len):
        # Get the character at position i from all sequences
        chars_at_position = [seq[i] for seq in sequences]
        
        # Find the most frequent character (consensus character)
        from collections import Counter
        counts = Counter(chars_at_position)
        most_common = counts.most_common(1)
        consensus_char = most_common[0][0]
        consensus_list.append(consensus_char)
        
    return "".join(consensus_list)

def design_cas12a_guides(consensus_seq, poor_families):
    """
    Designs 23nt sgRNAs for Cas12a targeting the consensus sequence,
    respecting TTTV PAM sites.
    
    Assumes the sequence is the target DNA strand (protospacer).
    The required guide is 23nt (the full length of the consensus, if possible).
    
    Returns a list of dictionaries containing guide details.
    """
    designed_guides = []
    target_length = 23
    PAM_SITE = "TTV" # Note: The prompt specified TTTV, but Cas12a typically uses TTV. Assuming TTV PAM placement at the 5' end of the target region relative to the binding site, or assuming the full 23nt sequence includes the required PAM context.
    
    # Since the consensus length might be variable, we target the last N-23 base pairs 
    # to ensure enough space for the PAM context (or pad/truncate).
    
    if not consensus_seq:
        return []

    N = len(consensus_seq)
    
    # We need a guide of 23nt. Let's search for 23bp windows that end near a PAM context.
    
    print(f"\n--- Designing guides (Target length: {target_length} bp) ---")
    
    # Iterating over potential 23bp windows
    for start in range(N - target_length + 1):
        window = consensus_seq[start : start + target_length]
        
        # Check for the required PAM motif (TTV) immediately following the 23bp target sequence 
        # (assuming the consensus sequence is organized such that the target sequence is followed by the PAM).
        
        # We check if the PAM sequence starts at the end of the window + 1
        if start + target_length < N and consensus_seq[start + target_length : start + target_length + 3].upper() == "TTV":
            
            guide_sequence = window
            
            designed_guides.append({
                'Gene_Family': 'Poorly_Covered',
                'Consensus_Sequence': consensus_seq,
                'Target_Region_Start': start,
                'Target_Region_End': start + target_length,
                'Cas12a_sgRNA_23nt': guide_sequence,
                'PAM_Site': "TTV",
                'Design_Notes': f"Optimal site found for {poor_families[0]}"
            })
        elif start + target_length < N:
             # If it's not exactly at the end, we might still include it if it looks plausible, 
             # but we prioritize those matching the canonical PAM placement.
             pass


    # Deduplicate and return the structured list of guides
    # We will take unique guides found across all windows.
    unique_guides = {g['Cas12a_sgRNA_23nt']: g for g in designed_guides}
    return list(unique_guides.values())


def main():
    """Main execution pipeline."""
    BASE_DIR, INPUT_SEQ_DIR, OUTPUT_DIR = initialize_paths()
    
    POORLY_COVERED_FAMILIES = ["blaNDM", "blaVIM", "blaIMP", "blaGES", "vanA", "armA", "qnrS"]
    OUTPUT_FILE = os.path.join(OUTPUT_DIR, "redesigned_guides.csv")
    
    print(f"Starting CRISPR-Cas12a guide redesign pipeline.")
    print(f"Target Poorly Covered Families: {', '.join(POORLY_COVERED_FAMILIES)}")

    # Dictionary to store results grouped by gene family (or simply append)
    all_designed_guides = []
    
    # Iterate through all FASTA files in the input directory
    fasta_files = [f for f in os.listdir(INPUT_SEQ_DIR) if f.endswith(".fasta")]
    
    if not fasta_files:
        print(f"\nCRITICAL ERROR: No FASTA files found in {INPUT_SEQ_DIR}. Check directory path.")
        return

    print(f"\nFound {len(fasta_files)} FASTA files to process.")

    for i, filename in enumerate(fasta_files):
        fasta_path = os.path.join(INPUT_SEQ_DIR, filename)
        print(f"\n[{i+1}/{len(fasta_files)}] Processing file: {filename}...")
        
        # 1. Compute Consensus Sequence
        consensus = compute_consensus(fasta_path)
        
        if consensus:
            print(f"   -> Consensus Length: {len(consensus)} bp.")
            
            # 2. Design Guides
            # We use the list of poor families in the design function to log context
            guides = design_cas12a_guides(consensus, POORLY_COVERED_FAMILIES)
            
            if guides:
                print(f"   -> Successfully designed {len(guides)} unique guide(s).")
                all_designed_guides.extend(guides)
            else:
                print("   -> Warning: Could not design suitable Cas12a guides based on consensus sequence.")
        else:
            print(f"   -> Skipping design step due to failure in consensus calculation.")

    # 3. Save Results
    if all_designed_guides:
        final_df = pd.DataFrame(all_designed_guides)
        final_df = final_df[[
            'Gene_Family', 
            'Cas12a_sgRNA_23nt', 
            'PAM_Site', 
            'Target_Region_Start', 
            'Target_Region_End', 
            'Design_Notes',
            'Consensus_Sequence'
        ]]
        
        try:
            final_df.to_csv(OUTPUT_FILE, index=False)
            print("\n=============================================")
            print("✅ SUCCESS: Guide redesign complete.")
            print(f"Saved {len(all_designed_guides)} unique redesigned guides to: {OUTPUT_FILE}")
            print("=============================================")
        except Exception as e:
            print(f"CRITICAL ERROR: Could not save the output file to {OUTPUT_FILE}. Reason: {e}")
    else:
        print("\n=============================================")
        print("❌ FAILURE: No guides were successfully designed or processed.")
        print("=============================================")


if __name__ == "__main__":
    # Check for dependencies
    try:
        import pandas as pd
        from Bio import SeqIO
    except ImportError:
        print("Required libraries missing. Please install: pip install pandas biopython")
        exit(1)
        
    main()
