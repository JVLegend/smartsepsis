import pandas as pd
import os
import sys
import io
import base64
from typing import List, Dict, Tuple

# --- Configuration and Paths ---
WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
OUTPUT_DIR = WORKING_DIR + "autoresearch_results/"

GUIDES_PATH = WORKING_DIR + "guides/*.tsv"
CONSERVATION_PATH = WORKING_DIR + "reports/conservation_analysis.tsv"
OUTPUT_FILE = OUTPUT_DIR + "multiplex_panel.csv"

def setup_directories():
    """Ensures the output directory exists."""
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        print(f"Created output directory: {OUTPUT_DIR}")

def load_data(file_path: str) -> pd.DataFrame:
    """
    Loads data with robust error handling.
    """
    print(f"\n--- Attempting to load data from: {file_path} ---")
    try:
        # Strip the glob pattern if it was used for the input
        if file_path.endswith("*.tsv"):
            base_path = os.path.join(WORKING_DIR, 'guides')
            # For guide loading, we assume iterating through files in the 'guides' subdirectory
            if "guides" in file_path and not os.path.isdir(base_path):
                 print(f"Error: Directory {base_path} not found.")
                 return pd.DataFrame()
            
            # Find all TSV files in the guides subdirectory
            guide_files = [os.path.join(base_path, f) for f in os.listdir(base_path) if f.endswith(".tsv")]
            if not guide_files:
                print(f"Error: No TSV files found in {base_path}")
                return pd.DataFrame()
            
            all_guide_data = []
            for i, f_path in enumerate(guide_files):
                print(f"Loading guides file {i+1}/{len(guide_files)}: {os.path.basename(f_path)}")
                try:
                    df = pd.read_csv(f_path, sep='\t')
                    all_guide_data.append(df)
                except Exception as e:
                    print(f"Warning: Failed to read {f_path}. Error: {e}")
            
            if all_guide_data:
                return pd.concat(all_guide_data, ignore_index=True)
            else:
                return pd.DataFrame()

        elif "conservation_analysis" in file_path:
            # Standard single file load
            return pd.read_csv(file_path, sep='\t')
        else:
            print(f"Error: Unknown file path format or type: {file_path}")
            return pd.DataFrame()
    except FileNotFoundError:
        print(f"ERROR: Required file not found at {file_path}")
        return pd.DataFrame()
    except pd.errors.EmptyDataError:
        print(f"ERROR: File is empty: {file_path}")
        return pd.DataFrame()
    except Exception as e:
        print(f"An unexpected error occurred loading {file_path}: {e}")
        return pd.DataFrame()

def optimize_multiplex_panel(guides_df: pd.DataFrame, conservation_df: pd.DataFrame) -> pd.DataFrame:
    """
    Performs greedy optimization to select the minimal set of guides
    that cover the maximum number of unique gene families.
    """
    print("\n--- Starting Multiplex Optimization (Set Cover Approximation) ---")
    
    if guides_df.empty or conservation_df.empty:
        print("Optimization aborted due to missing input data.")
        return pd.DataFrame()

    # 1. Determine the set of target families (Union of families present in both reports)
    guide_families = set(guides_df['gene_family'].unique())
    conservation_families = set(conservation_df['family'].unique())
    
    target_families = guide_families.union(conservation_families)
    print(f"Total unique target families identified: {len(target_families)}")
    
    if not target_families:
        print("No target families found in either dataset.")
        return pd.DataFrame()

    # 2. Greedy Set Cover Algorithm
    selected_guide_ids = set()
    covered_families = set()
    final_guides = []
    
    # Loop until all target families are covered, or we run out of useful guides
    while covered_families != target_families:
        best_guide_id = None
        max_new_coverage = 0
        
        # Iterate through all guides not yet selected
        potential_guides = guides_df[~guides_df['guide_id'].isin(selected_guide_ids)]
        
        if potential_guides.empty:
            print("Warning: All guides have been selected or no more unique guides available.")
            break

        # Find the guide that covers the maximum number of UNCOVERED families
        for index, row in potential_guides.iterrows():
            guide_id = row['guide_id']
            current_families = {row['gene_family']} # Assuming guide_family is the primary family link
            
            # Calculate the new coverage count
            new_coverage = len(current_families - covered_families)
            
            if new_coverage > max_new_coverage:
                max_new_coverage = new_coverage
                best_guide_id = guide_id
                best_guide_row = row
        
        # If no new coverage can be achieved, stop
        if best_guide_id is None or max_new_coverage == 0:
            break
            
        # Select the best guide
        selected_guide_ids.add(best_guide_id)
        final_guides.append(best_guide_row)
        
        # Update covered families
        new_families = {best_guide_row['gene_family']}
        covered_families.update(new_families)
        
        print(f"\rProgress: Selected {len(final_guides)} guides. Covered {len(covered_families)}/{len(target_families)} families.")

    print("\nOptimization Complete.")

    # 3. Compile the final output DataFrame
    if not final_guides:
        print("Could not select any multiplex panel.")
        return pd.DataFrame()
        
    multiplex_df = pd.DataFrame(final_guides)
    
    # Select and reorder key columns for the final output
    output_cols = ['guide_id', 'sequence', 'gene_family', 'pam_site']
    return multiplex_df[output_cols].drop_duplicates()


def main():
    """Main execution pipeline function."""
    print("="*80)
    print("CRISPR-Cas12a AMR Diagnostic Pipeline: Multiplex Optimization")
    print(f"Working Directory: {WORKING_DIR}")
    print("="*80)

    setup_directories()
    
    # 1. Load Guides Data
    guides_df = load_data(GUIDES_PATH)
    if guides_df.empty:
        print("\n!!! FATAL: Could not load guides data. Exiting.")
        return

    # 2. Load Conservation Analysis Data
    conservation_df = load_data(CONSERVATION_PATH)
    if conservation_df.empty:
        print("\n!!! WARNING: Conservation analysis data is missing or empty. Proceeding with guide-only optimization.")
    
    # 3. Run Optimization
    multiplex_df = optimize_multiplex_panel(guides_df, conservation_df)
    
    # 4. Save Output
    if not multiplex_df.empty:
        try:
            multiplex_df.to_csv(OUTPUT_FILE, index=False)
            print("\n" + "="*80)
            print(f"SUCCESS: Multiplex panel generated successfully.")
            print(f"Output saved to: {OUTPUT_FILE}")
            print("="*80)
        except Exception as e:
            print(f"\nERROR: Failed to write output file {OUTPUT_FILE}. Reason: {e}")
    else:
        print("\n!!! WARNING: No viable multiplex panel generated. Check input data.")

if __name__ == "__main__":
    main()
