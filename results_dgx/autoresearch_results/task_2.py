import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import base64
import io
import os
from typing import Dict, List

# Define required paths
WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
OUTPUT_DIR = os.path.join(WORKING_DIR, "autoresearch_results/")
INPUT_FILE = os.path.join(WORKING_DIR, "reports/evo2_scoring/functional_scores.tsv")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# TASK: B - Embedding comparison: read reports/evo2_scoring/functional_scores.tsv, compute pairwise distance matrix between all variants, save as CSV and generate an HTML heatmap using base64-encoded matplotlib.

def compute_pairwise_distance_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes the pairwise distance matrix based on the 'functional_distance' column.
    """
    print("--- Starting pairwise distance matrix computation ---")
    
    # 1. Extract the core distance data (functional_distance)
    distance_series = df['functional_distance'].values
    
    # 2. Calculate the pairwise distance matrix using scipy/numpy (or a pandas implementation)
    # Since we only rely on numpy/pandas, we use a direct pairwise calculation loop or a specialized library function.
    # Using numpy for efficiency
    N = len(distance_series)
    distance_matrix = np.zeros((N, N))
    
    for i in range(N):
        for j in range(i, N):
            distance_matrix[i, j] = distance_series[i]
            distance_matrix[j, i] = distance_series[i] # Assuming the input functional_distance represents the distance from a conceptual baseline/reference,
                                                        # but for pairwise variance calculation, we must assume the functional_distance column
                                                        # contains a property that can be used as a coordinate or dimension.
                                                        
            # Correction: If 'functional_distance' is a property value, and we need a pairwise distance, 
            # we must assume the *row index* defines the variant identity, and the distance calculation 
            # should be based on the actual distances between variants, not just the functional_distance column.
            # However, given the column name 'functional_distance' and the context of "Compute pairwise distance matrix,"
            # we assume this column itself holds a quantitative feature, and the distance between two variants 
            # (V1, V2) is simply the absolute difference of their functional_distance: |FD_V1 - FD_V2|.
            
            # Let's calculate the absolute difference pairwise.
            if i == j:
                distance_matrix[i, j] = 0.0 # Distance of a variant to itself is zero
            else:
                distance_matrix[i, j] = np.abs(distance_series[i] - distance_series[j])

    return pd.DataFrame(distance_matrix, index=df['variant'], columns=df['variant'])

def generate_heatmap_and_html(df_matrix: pd.DataFrame):
    """
    Generates the heatmap visualization and saves the results to an embedded HTML file.
    """
    print("--- Generating Heatmap and HTML Report ---")
    
    plt.figure(figsize=(12, 10))
    
    # Create the heatmap using seaborn
    # Using 'annot=False' for large matrices to save space, but showing structure is key.
    heatmap = sns.heatmap(
        df_matrix, 
        cmap="viridis", 
        annot=False, 
        cbar_kws={'label': 'Pairwise Functional Distance'}
    )
    plt.title("Pairwise Functional Distance Matrix of Variants (CRISPR-Cas12a)")
    plt.xlabel("Variant")
    plt.ylabel("Variant")
    
    # ---------------------------------------------------
    # Base64 Encoding
    # ---------------------------------------------------
    
    # Save the figure to a BytesIO object
    buf = io.BytesIO()
    plt.savefig(buf, format="png")
    plt.close() # Crucial to free memory
    
    # Encode the bytes to base64 string
    image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    
    # ---------------------------------------------------
    # Generating HTML Report
    # ---------------------------------------------------
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Pairwise Functional Distance Heatmap</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f4f4f9; }}
            .container {{ max-width: 1000px; margin: auto; padding: 20px; background: white; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
            h1 {{ color: #333; text-align: center; border-bottom: 2px solid #eee; padding-bottom: 10px; }}
            img {{ max-width: 100%; height: auto; display: block; margin: 20px auto; border: 1px solid #ddd; border-radius: 4px; }}
            p {{ color: #555; line-height: 1.6; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>CRISPR-Cas12a Variant Pairwise Functional Distance Analysis</h1>
            <p>This heatmap visualizes the pairwise absolute difference in the functional distance score between all analyzed resistance variants. Closer values (darker color, depending on colormap) indicate variants with similar functional effects relative to the baseline defined by the model.</p>
            <h2>Visualization</h2>
            <img src="data:image/png;base64,{image_base64}" alt="Pairwise Distance Heatmap">
            <p><strong>Methodology:</strong> The matrix is calculated as $|FD_{V_i} - FD_{V_j}|$, where $FD$ is the reported functional distance for variant $V$.</p>
            <p><strong>Summary:</strong> The resulting distance matrix quantifies the geometric separation of resistance mechanisms in the functional space.</p>
        </div>
    </body>
    </html>
    """
    
    output_html_path = os.path.join(OUTPUT_DIR, "pairwise_distance_heatmap.html")
    with open(output_html_path, "w") as f:
        f.write(html_content)
        
    print(f"\n[SUCCESS] Heatmap and HTML report saved successfully.")
    print(f"[PATH] Output HTML file: {output_html_path}")


def main():
    """
    Main execution function for Task B.
    """
    print("=========================================================")
    print("Task B: Pairwise Embedding Comparison Pipeline Start")
    print("=========================================================")

    try:
        # 1. Load the functional scores data
        print(f"Attempting to load data from: {INPUT_FILE}")
        df = pd.read_csv(INPUT_FILE, sep='\t')
        print(f"[INFO] Data loaded successfully. Number of records (variants): {len(df)}")

        if 'functional_distance' not in df.columns:
            raise ValueError("The required column 'functional_distance' is missing from the input file.")

        # 2. Compute the pairwise distance matrix
        distance_df = compute_pairwise_distance_matrix(df)
        
        # 3. Save the matrix to CSV
        csv_output_path = os.path.join(OUTPUT_DIR, "pairwise_distance_matrix.csv")
        distance_df.to_csv(csv_output_path, index=True)
        print(f"[SUCCESS] Pairwise distance matrix saved to: {csv_output_path}")

        # 4. Generate Heatmap and HTML Report
        generate_heatmap_and_html(distance_df)

    except FileNotFoundError:
        print("\n[ERROR] Failed to find the required input file.")
        print(f"[PATH] Expected file: {INPUT_FILE}")
        print("Please ensure the data path and filename are correct.")
    except ValueError as e:
        print(f"\n[ERROR] Data validation failed: {e}")
    except Exception as e:
        print(f"\n[CRITICAL ERROR] An unexpected error occurred during processing: {e}")
    finally:
        print("\n=========================================================")
        print("Task B: Pairwise Embedding Comparison Pipeline Finished.")
        print("=========================================================")

if __name__ == "__main__":
    main()
