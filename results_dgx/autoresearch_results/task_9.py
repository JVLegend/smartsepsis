import pandas as pd
import os
import numpy as np
import io
import base64
from typing import List, Dict

# Define required absolute paths
BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
GUIDES_PATH = os.path.join(BASE_DIR, "guides")
OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "guide_quality.csv")

# --- Helper Functions ---

def calculate_gc_content(sequence: str) -> float:
    """Calculates the GC percentage of a given DNA sequence."""
    if not sequence:
        return 0.0
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    gc_content = (g_count + c_count) / len(sequence)
    return gc_content * 100

def calculate_tm_nearest_neighbor(sequence: str) -> float:
    """
    Calculates the theoretical Melting Temperature (Tm) using a simplified
    nearest-neighbor model approximation (Placeholder for complex calculation).
    
    Note: A true nearest-neighbor calculation requires specific empirical
    tables for optimal accuracy. This function simulates the calculation 
    process using a standard simplified formula structure.
    """
    if len(sequence) < 10:
        return 30.0  # Too short for accurate estimation
    
    # Example simplification: Use a general formula that correlates length and GC content
    gc_percent = calculate_gc_content(sequence)
    
    # Basic estimate: Tm = 2*(A+T) + 4*(G+C) (For oligonucleotides <= 20 bases)
    # If the guide is longer, a linear scaling factor (L/2) must be applied.
    L = len(sequence)
    
    # Let's use a weighted average approach simulating the nearest neighbor effect
    tm = (2 * ((100 - gc_percent) / 4) + (4 * gc_percent / 4)) + (L / 2.0)
    
    # Ensure a minimum reasonable Tm
    return max(30.0, tm)

def estimate_hairpin_potential(sequence: str) -> float:
    """
    Estimates the minimum free energy (or potential) for hairpin formation.
    Returns a relative score (negative is better/stronger potential).
    """
    # In a real scenario, this would use algorithms like mfold or ViennaRNA.
    # For demonstration, we use a length-based and composition-based heuristic.
    
    # Penalty function: length and low diversity increase potential
    potential = -1.0 * (len(sequence) / 20.0) 
    potential -= 0.05 * np.std([int(s) for s in sequence])
    
    return potential

# --- Main Processing Function ---

def process_guide_quality_analysis():
    """
    Reads guide sequences, computes GC content, Tm, and hairpin potential,
    and flags guides based on defined quality criteria.
    """
    print("=" * 50)
    print("CRISPR-Cas12a Guide Quality Analysis Pipeline Started.")
    print(f"Input guides directory: {GUIDES_PATH}")
    print(f"Output file: {OUTPUT_FILE}")
    print("=" * 50)
    
    all_data_records: List[Dict] = []
    
    try:
        # Check if the guides directory exists
        if not os.path.exists(GUIDES_PATH):
            raise FileNotFoundError(f"The guides directory was not found at: {GUIDES_PATH}")

        guide_files = [f for f in os.listdir(GUIDES_PATH) if f.endswith('.tsv')]
        print(f"Found {len(guide_files)} guide TSV files to process.")

        for i, filename in enumerate(guide_files):
            file_path = os.path.join(GUIDES_PATH, filename)
            print(f"\n[PROGRESS] Processing file {i+1}/{len(guide_files)}: {filename}...")

            try:
                # Assuming the TSV structure matches the specified columns: 
                # guide_id, sequence, gene_family, pam_site
                df = pd.read_csv(file_path, sep='\t')
                
                if df.empty:
                    print(f"WARNING: {filename} was empty. Skipping.")
                    continue

                # Apply analysis functions row-wise
                df['GC_percent'] = df['sequence'].apply(calculate_gc_content)
                df['Tm_predicted'] = df['sequence'].apply(calculate_tm_nearest_neighbor)
                df['Hairpin_potential'] = df['sequence'].apply(estimate_hairpin_potential)
                
                # --- Flagging Logic ---
                
                # Criteria: GC<30% or >70%, Tm<50 or >70
                def flag_guide(row):
                    gc_flag = 'PASS'
                    if row['GC_percent'] < 30.0 or row['GC_percent'] > 70.0:
                        gc_flag = 'FAIL (GC)'
                        
                    tm_flag = 'PASS'
                    if row['Tm_predicted'] < 50.0 or row['Tm_predicted'] > 70.0:
                        tm_flag = 'FAIL (Tm)'
                        
                    return f"{gc_flag}; {tm_flag}"

                df['Quality_Flags'] = df.apply(flag_guide, axis=1)

                # Select and rename final columns
                results_df = df[[
                    'guide_id', 
                    'sequence', 
                    'gene_family', 
                    'GC_percent', 
                    'Tm_predicted', 
                    'Hairpin_potential', 
                    'Quality_Flags'
                ]].copy()
                
                results_df = results_df.rename(columns={
                    'GC_percent': 'GC_percent (%)',
                    'Tm_predicted': 'Tm_predicted (C)',
                    'Hairpin_potential': 'Hairpin_Potential_Score',
                })
                
                all_data_records.append(results_df)

            except pd.errors.EmptyDataError:
                print(f"ERROR: Could not read {filename}. File might be empty or corrupted.")
            except Exception as e:
                print(f"CRITICAL ERROR processing {filename}: {e}")
                
        # Concatenate all processed dataframes
        if all_data_records:
            final_results = pd.concat(all_data_records, ignore_index=True)
            
            # Save the final dataframe
            final_results.to_csv(OUTPUT_FILE, index=False)
            
            print("\n" + "=" * 50)
            print("TASK COMPLETE.")
            print(f"Successfully processed {len(all_data_records)} batch files.")
            print(f"Results saved to: {OUTPUT_FILE}")
            print("=" * 50)
        else:
            print("\n[FAILURE] No data was successfully processed. Check input files and paths.")

    except FileNotFoundError as e:
        print(f"\n[FATAL ERROR] Setup failed: {e}")
    except Exception as e:
        print(f"\n[FATAL ERROR] An unexpected error occurred: {e}")


if __name__ == "__main__":
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    process_guide_quality_analysis()