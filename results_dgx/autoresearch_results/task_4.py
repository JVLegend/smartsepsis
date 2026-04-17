import os
import glob
import csv
import itertools
from Bio import SeqIO
from Bio.Seq import Seq

# --- Configuration ---
BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
GUIDES_DIR = os.path.join(BASE_DIR, "guides/")
SEQUENCES_DIR = os.path.join(BASE_DIR, "sequences/")
OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results/")
CROSS_REACTIVITY_FILE = os.path.join(OUTPUT_DIR, "cross_reactivity.csv")
TARGET_LENGTH = 23
MAX_MISMATCHES = 3

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# TASK: Cross-reactivity analysis: read all guide TSV files from guides/ dir. 
# For each guide, check if its 23nt sequence appears (with <=3 mismatches) 
# in FASTA files of OTHER gene families. Report off-target hits as cross_reactivity.csv.

def parse_guides(guides_dir):
    """Reads and structures guide data from TSV files."""
    guide_data = []
    print(f"Reading guide files from: {guides_dir}")
    
    guide_paths = glob.glob(os.path.join(guides_dir, "*.tsv"))
    if not guide_paths:
        print("Warning: No guide TSV files found. Exiting.")
        return None

    for file_path in guide_paths:
        try:
            print(f"Processing guide file: {file_path}")
            with open(file_path, 'r') as f:
                reader = csv.reader(f)
                # Skip header
                next(reader)
                for row in reader:
                    if len(row) >= 4:
                        guide_data.append({
                            'guide_id': row[0],
                            'sequence': row[1].strip().upper(),
                            'gene_family': row[2].strip().upper(),
                            'pam_site': row[3].strip().upper()
                        })
        except Exception as e:
            print(f"Error reading or processing guide file {file_path}: {e}")
    
    print(f"Successfully loaded {len(guide_data)} guides.")
    return guide_data

def load_fasta_sequences(sequences_dir):
    """Loads all FASTA files and stores sequences mapped by gene family."""
    print(f"Loading all FASTA sequences from: {sequences_dir}")
    
    all_sequences = {}
    fasta_paths = glob.glob(os.path.join(sequences_dir, "*.fasta"))
    
    for fasta_file in fasta_paths:
        try:
            # Use BioPython SeqIO for reliable FASTA parsing
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq = str(record.seq).upper()
                accession = record.id
                # Assuming the gene family is derivable from the filename or metadata
                # For this simulation, we'll use a simple heuristic: the directory structure might encode the family.
                # We'll use the filename stem as a proxy for 'gene_family' if no clear metadata is available.
                # NOTE: In a real pipeline, the family would be explicitly provided.
                
                # Extract a pseudo-family based on the parent directory or file prefix
                pseudo_family = os.path.basename(os.path.dirname(fasta_file))
                
                if pseudo_family not in all_sequences:
                    all_sequences[pseudo_family] = []
                
                all_sequences[pseudo_family].append({
                    'accession': accession,
                    'sequence': seq,
                    'family': pseudo_family # Using the directory name as the family identifier
                })
        except Exception as e:
            print(f"Error loading FASTA file {fasta_file}: {e}")
    
    return all_sequences

def hamming_distance(seq1, seq2):
    """Calculates the number of mismatches (Hamming distance). Assumes equal length."""
    if len(seq1) != len(seq2):
        # Should not happen if slicing is done correctly, but good practice.
        return -1 
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def search_cross_reactivity(guide_data, sequence_map):
    """
    Performs the core cross-reactivity search.
    Checks if the guide sequence matches any target sequence from a different family 
    with <= 3 mismatches.
    """
    print("\n--- Starting Cross-Reactivity Search ---")
    hits = []
    total_guides = len(guide_data)
    search_counter = 0

    # Iterate through every guide
    for guide in guide_data:
        guide_seq = guide['sequence']
        guide_family = guide['gene_family']
        guide_id = guide['guide_id']

        cross_reactive_hits = []
        
        # Iterate through every gene family in the loaded sequences
        for target_family, targets in sequence_map.items():
            
            # Rule: Exclude searches within the same gene family
            if target_family == guide_family:
                continue
            
            # Iterate through every sequence within the target family
            for target_record in targets:
                target_seq = target_record['sequence']
                target_accession = target_record['accession']

                # Use a sliding window approach to check all possible starting positions
                # The search must maintain the guide length (23nt)
                for i in range(len(target_seq) - TARGET_LENGTH + 1):
                    target_window = target_seq[i:i + TARGET_LENGTH]
                    
                    # Calculate mismatch count
                    mismatches = hamming_distance(guide_seq, target_window)
                    
                    if mismatches <= MAX_MISMATCHES:
                        hit = {
                            'guide_id': guide_id,
                            'guide_sequence': guide_seq,
                            'guide_family': guide_family,
                            'hit_accession': target_accession,
                            'hit_family': target_family,
                            'hit_sequence': target_window,
                            'mismatches': mismatches,
                            'hit_position': i # Starting position in the target sequence
                        }
                        
                        # Avoid duplicate reporting for the same guide/target pair
                        if hit not in cross_reactive_hits:
                            cross_reactive_hits.append(hit)
                
            search_counter += 1
            
        if cross_reactive_hits:
            hits.extend(cross_reactive_hits)
        
        print(f"Processed Guide {guide_id}: Found {len(cross_reactive_hits)} hits.")
        
        # Progress tracking
        if (search_counter % 10 == 0) or (search_counter == 1):
             print(f"--- Progress: {search_counter}/{len(sequence_map)} target families checked ---")


    return hits

def save_results(hits, output_path):
    """Saves the list of hits to the specified CSV file."""
    if not hits:
        print("\nNo cross-reactivity hits found to report.")
        return

    fieldnames = [
        'guide_id', 
        'guide_sequence', 
        'guide_family', 
        'hit_accession', 
        'hit_family', 
        'hit_sequence', 
        'mismatches', 
        'hit_position'
    ]
    
    print(f"\nSaving {len(hits)} cross-reactivity hits to: {output_path}")
    try:
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(hits)
        print("Cross-reactivity analysis complete. Results saved successfully.")
    except IOError as e:
        print(f"FATAL ERROR writing output file {output_path}: {e}")


# --- Main Execution ---
def main():
    """Orchestrates the entire diagnostic pipeline."""
    print("==============================================================")
    print("  CRISPR-Cas12a Cross-Reactivity Diagnostic Pipeline START ")
    print("==============================================================")
    
    # 1. Load Guides
    guide_data = parse_guides(GUIDES_DIR)
    if guide_data is None:
        return
    
    # 2. Load Target Sequences (FASTA)
    sequence_map = load_fasta_sequences(SEQUENCES_DIR)
    if not sequence_map:
        print("No target FASTA sequences found. Aborting analysis.")
        return
        
    # 3. Run Cross-Reactivity Search
    cross_reactive_hits = search_cross_reactivity(guide_data, sequence_map)
    
    # 4. Save Results
    save_results(cross_reactive_hits, CROSS_REACTIVITY_FILE)

if __name__ == "__main__":
    main()
