import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import os

def compute_and_plot_distance_matrix(input_path: str, output_csv: str, output_html: str):
    """
    Reads functional pairwise scores, computes the distance matrix, saves it to CSV,
    and generates a base64-encoded HTML heatmap.
    """
    print(f"--- Starting Embedding Comparison Workflow ---")

    # 1. Load Data
    try:
        df = pd.read_csv(input_path, sep='\t')
        print(f'Successfully loaded data from: {input_path}')
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        return
    except Exception as e:
        print(f"An error occurred during file loading: {e}")
        return

    print(f'Columns: {df.columns.tolist()}')

    # Ensure 'functional_distance' column exists and is numeric
    if 'functional_distance' not in df.columns:
        print("Error: 'functional_distance' column required but not found.")
        return

    # 2. Prepare the Distance Matrix Data
    # We need the matrix of distances for all pairs of unique variants.
    
    # Filter to only include the relevant columns for matrix calculation
    # We assume (variant, reference) pairs define the functional distance.
    df_scores = df[['variant', 'reference', 'functional_distance']].copy()

    # Identify all unique variants involved
    all_variants = pd.unique(df_scores[['variant', 'reference']].values.ravel('K'))
    
    # Create a dictionary mapping (v1, v2) -> distance
    distance_data = {}
    
    for index, row in df_scores.iterrows():
        v1 = row['variant']
        v2 = row['reference']
        distance = row['functional_distance']
        
        # Standardize key order for symmetry: ensure (A, B) is treated the same as (B, A)
        pair = tuple(sorted((v1, v2)))
        
        # Use a consistent key format (e.g., v1_v2 where v1 <= v2)
        key = f"{pair[0]}_{pair[1]}"
        
        # Since the input might contain both (A, B) and (B, A) rows, we take the mean/last encountered value.
        # We'll treat the distance as given in the row, assuming the input is structured correctly.
        distance_data[(v1, v2)] = distance

    # 3. Compute the Symmetric Distance Matrix
    
    # Initialize the matrix
    matrix_size = len(all_variants)
    dist_matrix = np.full((matrix_size, matrix_size), np.nan)
    
    # Map variant names to indices
    variant_to_index = {v: i for i, v in enumerate(all_variants)}

    # Populate the matrix
    for (v1, v2), dist in distance_data.items():
        i1 = variant_to_index[v1]
        i2 = variant_to_index[v2]
        
        # Populate both directions to ensure symmetry (D[i, j] = D[j, i])
        dist_matrix[i1, i2] = dist
        dist_matrix[i2, i1] = dist
    
    # Set diagonal to 0 (distance from a variant to itself)
    np.fill_diagonal(dist_matrix, 0.0)

    # 4. Save Matrix to CSV
    df_matrix = pd.DataFrame(dist_matrix, index=all_variants, columns=all_variants)
    df_matrix.index.name = 'Variant'
    df_matrix.columns.name = 'Variant'
    
    df_matrix.to_csv(output_csv)
    print(f'\n[SUCCESS] Distance matrix saved to {output_csv} (CSV format).')

    # 5. Generate Heatmap
    plt.figure(figsize=(15, 13))
    
    # Use correlation plot style since we are plotting distances (which are inherently non-correlated values)
    # Use a divergent colormap for better visualization of 0 (minimum distance)
    sns.heatmap(
        df_matrix, 
        annot=False, # Annotation might clutter a large matrix
        cmap="viridis", 
        linewidths=.5, 
        cbar_kws={'label': 'Functional Distance'}
    )
    plt.title('Pairwise Functional Distance Heatmap of Variants', fontsize=16)
    plt.xlabel('Variant (Pair)')
    plt.ylabel('Variant (Pair)')
    plt.tight_layout()
    
    # 6. Encode and Save HTML
    buf = BytesIO()
    plt.savefig(buf, format='png')
    plt.close()
    
    # Base64 encoding
    with open("heatmap.png", "rb") as f:
        img_base64 = base64.b64encode(f.read()).decode('utf-8')
    
    # Create a simple HTML structure incorporating the base64 image
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Functional Distance Heatmap</title>
        <style>
            body {{ font-family: Arial, sans-serif; text-align: center; padding: 20px; }}
            img {{ max-width: 90%; height: auto; border: 1px solid #ccc; }}
        </style>
    </head>
    <body>
        <h1>Pairwise Functional Distance Heatmap</h1>
        <p>This heatmap visualizes the functional distance between all identified variants.</p>
        <img src="data:image/png;base64,{img_base64}" alt="Functional Distance Heatmap">
    </body>
    </html>
    """
    
    with open(output_html, "w") as f:
        f.write(html_content)
        
    # Cleanup temporary image file
    os.remove("heatmap.png")
    
    print(f'\n[SUCCESS] Heatmap saved as base64-encoded HTML to {output_html}.')
    print("--- Workflow Complete ---")


if __name__ == "__main__":
    INPUT_FILE = "reports/evo2_scoring/functional_scores.tsv"
    OUTPUT_CSV = "functional_distance_matrix.csv"
    OUTPUT_HTML = "functional_distance_heatmap.html"
    
    compute_and_plot_distance_matrix(INPUT_FILE, OUTPUT_CSV, OUTPUT_HTML)
