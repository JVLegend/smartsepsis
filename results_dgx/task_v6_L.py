import pandas as pd
import glob
import os
import warnings
import numpy as np

# Suppress UserWarnings often encountered in bioinformatics workflows
warnings.filterwarnings('ignore', category=UserWarning)

def load_data(path_pattern, sep=None):
    """Loads data using pandas, handling TSV/CSV based on separator."""
    sep = sep if sep is not None else '\\t'
    try:
        df = pd.read_csv(path_pattern, sep=sep)
        print(f"Successfully loaded data from {path_pattern}")
        return df
    except FileNotFoundError:
        print(f"ERROR: File not found at {path_pattern}. Skipping.")
        return pd.DataFrame()
    except Exception as e:
        print(f"An error occurred while reading {path_pattern}: {e}")
        return pd.DataFrame()

def get_multiplex_panel(guide_dir, conservation_file):
    """
    Performs multiplex optimization: selects a minimal set of guides 
    covering maximum unique gene families with exact matches.
    """

    # --- 1. Load Data ---
    
    # 1.1 Load Conservation Data (Tsv)
    conservation_path = os.path.join(conservation_file)
    conservation_df = load_data(conservation_path)
    if conservation_df.empty:
        return pd.DataFrame()

    # 1.2 Load Guides (Tsv files)
    guide_files = glob.glob(os.path.join(guide_dir, '*.tsv'))
    guide_list = []
    for i, f_path in enumerate(guide_files):
        df = load_data(f_path)
        if not df.empty:
            df['source_file'] = os.path.basename(f_path)
            guide_list.append(df)
    
    all_guides_df = pd.concat(guide_list, ignore_index=True)
    if all_guides_df.empty:
        return pd.DataFrame()

    # --- 2. Filter Conservation Data ---
    # We only care about 'exact' matches
    exact_match_df = conservation_df[conservation_df['match'] == 'exact'].copy()
    
    # Clean up column names in the conserved data for easier joining
    if 'family' not in exact_match_df.columns or 'gene' not in exact_match_df.columns:
         # If the gene name is not explicitly present, we might need to rely on a proxy 
         # relationship (assuming 'variant' or 'accession' links to 'gene').
         # Given the schema provided, we assume that 'family' is the key metric 
         # and that the guiding mechanism relates the guide's 'gene' to the conserved 'family'.
         print("Warning: 'gene' column linking guides and conservation is critical. Assuming direct relationship between guide.gene and conserved context.")
         pass

    # --- 3. Multiplex Selection Logic (Greedy Algorithm) ---

    print("\n--- Starting Multiplex Selection ---")
    
    # 3.1 Identify unique required families
    required_families = exact_match_df['family'].unique()
    print(f"Total unique target families requiring coverage: {len(required_families)}")

    # 3.2 Create a mapping structure: family -> list of relevant guides
    # A guide is relevant if its 'gene' field matches or relates to the family/gene context.
    
    # Create a composite key: (guide gene, family name)
    guide_family_coverage = []
    for index, row in all_guides_df.iterrows():
        # Simple heuristic: Check if the guide's gene name appears in the conserved context.
        # This links the guide to the family it is intended to target.
        
        # Identify families related to the guide's gene
        related_families = exact_match_df[
            (exact_match_df['family'] == row['gene']) # Simple direct match
        ]['family'].unique().tolist()
        
        # If the family names don't match the gene name, try linking by accession or context (more complex, skipped for general robust script)
        if not related_families:
             # Fallback: If the guide gene isn't listed in the family column, skip for now.
             pass

        # Append the guide and all families it covers
        for family in related_families:
            guide_family_coverage.append({
                'guide_index': index,
                'guide_gene': row['gene'],
                'spacer_seq': row['spacer_seq'],
                'guide_identifier': row.get('rank', row.get('crRNA_AsCas12a', f"idx_{index}")),
                'covered_family': family
            })

    coverage_df = pd.DataFrame(guide_family_coverage)
    
    if coverage_df.empty:
        print("No guide-family coverage established. Outputting empty panel.")
        return pd.DataFrame()

    # 3.3 Greedy Selection Process
    
    selected_guide_indices = set()
    selected_guides = {} # {unique_guide_id: selected_row}
    covered_families = set()
    
    # While we still have uncovered families and available guides
    while covered_families != set(required_families):
        best_guide_info = None
        max_new_coverage = -1
        
        # Candidates are guides that haven't been selected yet
        candidate_coverage = coverage_df[~coverage_df['guide_index'].isin(list(selected_guide_indices))]
        
        if candidate_coverage.empty:
            break # No more guides to select
        
        for _, row in candidate_coverage.iterrows():
            current_guide_id = row['guide_identifier']
            current_guide_index = row['guide_index']
            
            # Calculate how many *new* families this guide covers
            newly_covered_families = set(row['covered_family']).difference(covered_families)
            new_coverage_count = len(newly_covered_families)
            
            # Selection Criteria: 
            # 1. Maximize new coverage.
            # 2. Tiebreaker: Prefer guides that are generally robust (e.g., high score, though score column is complex/missing).
            # We stick to the primary goal: max new coverage.
            
            if new_coverage_count > max_new_coverage:
                max_new_coverage = new_coverage_count
                best_guide_info = row
            elif new_coverage_count == max_new_coverage and best_guide_info is not None:
                 # Simple tiebreaker: keep the current choice or rely on index order
                 pass 

        if best_guide_info is None or max_new_coverage == 0:
            # No remaining guides can cover new families
            break

        # Select the best guide
        selected_guide_indices.add(best_guide_info['guide_index'])
        
        # Get the unique identifier for the selected guide
        guide_id = best_guide_info['guide_identifier']
        
        # Store the guide data (using the full row data)
        if guide_id not in selected_guides:
             # We must take the full, original row data
            original_guide_row = all_guides_df.loc[best_guide_info['guide_index']]
            selected_guides[guide_id] = original_guide_row

        # Update coverage
        new_families = set(best_guide_info['covered_family'])
        covered_families.update(new_families)

    # 3.4 Compile Final Panel
    final_df = pd.DataFrame(selected_guides.values())
    
    # Clean up the output for presentation (e.g., dropping internal index/source columns)
    # Assuming we want the core guide information: rank, gene, spacer_seq, etc.
    required_cols = [
        'rank', 'gene', 'spacer_seq', 'pam_seq', 'position', 'strand', 
        'spacer_length', 'score', 'gc', 'homopolymer', 'self_comp', 'poly_t', 
        'rel_position', 'crRNA_LbCas12a', 'crRNA_AsCas12a'
    ]
    # Filter to only columns that actually exist in the data
    final_columns = [col for col in required_cols if col in final_df.columns]
    
    multiplex_panel = final_df[final_columns].copy()
    
    print("\n--- Selection Complete ---")
    print(f"Target Families Covered: {covered_families.intersection(required_families)}")
    print(f"Total unique guides selected for multiplex panel: {len(multiplex_panel)}")
    
    return multiplex_panel


if __name__ == "__main__":
    # --- Configuration ---
    GUIDES_DIR = "guides"
    CONSERVATION_FILE = "reports/conservation_analysis.tsv"
    OUTPUT_FILE = "multiplex_panel.csv"

    # --- Execution ---
    panel = get_multiplex_panel(GUIDES_DIR, CONSERVATION_FILE)
    
    # --- Output ---
    if not panel.empty:
        # Output the final DataFrame to a CSV file
        panel.to_csv(OUTPUT_FILE, index=False)
        print(f"\n✅ Success! Multiplex panel saved to {OUTPUT_FILE}")
    else:
        print("\n❌ Failed to generate panel due to insufficient or malformed data.")
