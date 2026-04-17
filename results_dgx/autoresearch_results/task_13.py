import pandas as pd
import os
from collections import defaultdict

def build_resistance_network(root_dir: str):
    """
    Reads card_new_variants.csv and builds an adjacency matrix/edge list
    representing gene families connected by shared drug classes.

    Args:
        root_dir: The base directory path.
    """
    input_path = os.path.join(root_dir, "reports", "card_new_variants.csv")
    output_dir = os.path.join(root_dir, "autoresearch_results")
    output_file = os.path.join(output_dir, "network_adjacency.csv")

    print(f"# TASK: Resistance mechanism network: read card_new_variants.csv. Build adjacency matrix of gene families connected by shared drug classes.")

    # 1. Setup directories
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory ensured at: {output_dir}")
    except OSError as e:
        print(f"Error creating output directory {output_dir}: {e}")
        return

    # 2. Load Data
    try:
        df = pd.read_csv(input_path)
        print(f"Successfully loaded data from {input_path}. Shape: {df.shape}")
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}. Aborting.")
        return
    except Exception as e:
        print(f"An unexpected error occurred during file loading: {e}. Aborting.")
        return

    # 3. Preprocessing and Cleaning
    # Normalize drug_classes into lists of unique drugs
    df['drug_classes_list'] = df['drug_classes'].str.split(r'[,\s]+').apply(lambda x: [d.strip() for d in x if d.strip()]).tolist()

    print(f"Preprocessing complete. Unique families identified: {df['pipeline_family'].nunique()}")

    # Dictionary to hold connections: {Family A: {Family B: {Shared Drugs: set()}}}
    connections = defaultdict(lambda: defaultdict(set))

    # 4. Graph Building Logic
    families = df['pipeline_family'].unique()
    
    print("Building pairwise connections between gene families...")

    # Optimized approach: Group by family and then compare drug sets
    for i in range(len(families)):
        family_A = families[i]
        
        # Filter rows belonging to the current family A
        df_A = df[df['pipeline_family'] == family_A].copy()
        
        # For each pair (A, B) where B > A (to avoid duplicate edges)
        for j in range(i + 1, len(families)):
            family_B = families[j]
            
            # Filter rows belonging to family B
            df_B = df[df['pipeline_family'] == family_B].copy()

            # Check all combinations of rows (A, B) to find shared drugs
            # This logic is complex if we assume drug classes are row-dependent,
            # but since the task implies overall family connection, we aggregate all drug classes for the family.

            # Aggregate all unique drugs for Family A and Family B
            drugs_A = set()
            for drugs in df_A['drug_classes_list']:
                drugs_A.update(drugs)
            
            drugs_B = set()
            for drugs in df_B['drug_classes_list']:
                drugs_B.update(drugs)
            
            # Find the intersection (shared drug classes)
            shared_drugs = drugs_A.intersection(drugs_B)

            if shared_drugs:
                # Store the connection (undirected graph)
                for drug in shared_drugs:
                    connections[family_A][family_B].add(drug)
                connections[family_B][family_A].add(drug) # Ensure bidirectional connection


    # 5. Formatting Output
    output_records = []
    for family_A in connections:
        for family_B, drugs in connections[family_A].items():
            if family_A < family_B: # Ensure consistent ordering and only list A < B
                # Create a single record for all shared connections between A and B
                num_shared = len(drugs)
                shared_drugs_str = ', '.join(sorted(list(drugs)))
                output_records.append({
                    'Family_A': family_A,
                    'Family_B': family_B,
                    'Shared_Drug_Classes': shared_drugs_str,
                    'Connection_Count': num_shared
                })

    # 6. Saving Output
    df_output = pd.DataFrame(output_records)
    
    try:
        df_output.to_csv(output_file, index=False)
        print(f"\nSUCCESS: Resistance network adjacency matrix saved successfully to: {output_file}")
    except Exception as e:
        print(f"\nCRITICAL ERROR: Could not save output file to {output_file}. Reason: {e}")


# --- Main Execution ---
if __name__ == "__main__":
    # Define the root directory as required
    ROOT_WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
    build_resistance_network(ROOT_WORKING_DIR)
