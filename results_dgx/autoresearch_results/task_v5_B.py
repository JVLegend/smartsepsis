import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import base64
import io
import os

# --- Configuration ---
BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results/")

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define input file paths
FUNCTIONAL_SCORES_PATH = os.path.join(BASE_DIR, "reports/evo2_scoring/functional_scores.tsv")
DISTANCE_MATRIX_OUTPUT_PATH = os.path.join(OUTPUT_DIR, "variant_functional_distance_matrix.csv")
HEATMAP_HTML_OUTPUT_PATH = os.path.join(OUTPUT_DIR, "functional_distance_heatmap.html")

# TASK: B - Embedding comparison: read reports/evo2_scoring/functional_scores.tsv, compute pairwise distance matrix between all variants using functional_distance values, save as CSV and generate an HTML heatmap using base64-encoded matplotlib.

print("--- Starting Task B: Embedding Comparison ---")

# 1. Load functional scores data
try:
    print(f"Reading functional scores from: {FUNCTIONAL_SCORES_PATH}")
    # Reading the data, assuming tab separation as specified
    df = pd.read_csv(FUNCTIONAL_SCORES_PATH, sep='\t')
    print("Successfully loaded functional scores.")
    print(f"Columns: {df.columns.tolist()}")
except FileNotFoundError:
    print(f"Error: Functional scores file not found at {FUNCTIONAL_SCORES_PATH}. Aborting Task B.")
    exit()
except Exception as e:
    print(f"An unexpected error occurred while reading the functional scores: {e}. Aborting Task B.")
    exit()

# 2. Compute Pairwise Distance Matrix
try:
    # Identify unique variants (variants are the targets for the matrix)
    # We assume the variants are defined by the 'variant' column in this context.
    variants = df['variant'].unique()
    print(f"\nFound {len(variants)} unique variants to analyze.")

    # Initialize an empty DataFrame to hold the matrix data
    distance_data = {v: [] for v in variants}

    # Use numpy/pandas functionality to compute the distance matrix
    # We must use functional_distance for this calculation.
    
    # Create a pivot table structure where rows are variants and columns are variants,
    # filled with their functional_distance.
    distance_df = df.pivot(index='variant', columns='reference', values='functional_distance')
    
    # Fill NaN values (which occur when a pair is missing in the input data) with a placeholder 
    # or handle them according to assumptions. For a true distance matrix, we must ensure all pairs exist.
    # A common simplification is to calculate the average or use a placeholder if the pair is missing.
    # Given the prompt implies a complete matrix calculation, we will fill NaNs with NaN 
    # and calculate the full matrix ensuring symmetry.
    
    # For simplicity and robustness, we enforce the pairing using the list of variants.
    distance_matrix = pd.DataFrame(index=variants, columns=variants)
    
    for i in range(len(variants)):
        for j in range(len(variants)):
            var1 = variants[i]
            var2 = variants[j]
            
            # Try to extract the distance from the original dataframe
            # The distance is stored when var1 is the 'variant' and var2 is the 'reference',
            # OR vice versa.
            
            # Check the existing pairs in the DataFrame
            match = df[
                ((df['variant'] == var1) & (df['reference'] == var2)) |
                ((df['variant'] == var2) & (df['reference'] == var1))
            ]
            
            if not match.empty:
                # We take the mean distance if multiple entries exist for the pair
                distance = match['functional_distance'].mean()
            else:
                # If the pair does not exist in the dataset, assign NaN or 0 depending on context. 
                # Using np.nan for robustness.
                distance = np.nan
            
            distance_matrix.iloc[i, j] = distance

    # Save the resulting distance matrix
    distance_matrix.to_csv(DISTANCE_MATRIX_OUTPUT_PATH)
    print(f"\nSuccessfully computed and saved pairwise distance matrix to: {DISTANCE_MATRIX_OUTPUT_PATH}")

except KeyError as e:
    print(f"Error: Missing required column for distance computation: {e}. Check the headers of functional_scores.tsv.")
except Exception as e:
    print(f"An error occurred during matrix computation: {e}")

# 3. Generate Heatmap Visualization
try:
    print("\nGenerating heatmap visualization...")
    
    # Set up the matplotlib figure
    plt.figure(figsize=(20, 18))
    
    # Use seaborn to generate the heatmap
    # We use a symmetric representation, so we trim the diagonal warning/plot
    mask = np.triu(np.ones_like(distance_matrix, dtype=bool))
    
    # Only plot the upper triangular part and set lower values to NaN to make it look clean
    # We must ensure the diagonal is also handled (distances to self should be 0 or NaN)
    
    # Create a copy for plotting
    plot_matrix = distance_matrix.copy()
    # Set the upper triangle and diagonal to NaN so they are not colored (optional, but cleaner)
    plot_matrix = plot_matrix.where(~mask)
    
    # Calculate the maximum distance for consistent coloring
    vmax = np.max(distance_matrix.dropna())
    
    # Use a diverging color map (like 'viridis_r' or 'cividis')
    # We use 'viridis' as a standard choice for distance heatmaps.
    sns.heatmap(
        plot_matrix, 
        annot=False, # Annotation takes too much space for 42 variants
        fmt=".1f", 
        cmap="viridis", 
        cbar_kws={'label': 'Functional Distance'}, 
        mask=mask, # Masking is preferred way for upper triangle plotting
        linewidths=.5, 
        linecolor='lightgray'
    )
    
    plt.title("Pairwise Functional Distance Matrix of Gene Variants (CRISPR-Cas12a)", fontsize=16)
    plt.xlabel("Variant", fontsize=14)
    plt.ylabel("Variant", fontsize=14)
    plt.xticks(rotation=90, ha='right')
    plt.yticks(rotation=0)
    
    # Save the figure to an in-memory buffer (BytesIO)
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close()
    
    # Encode the image buffer to base64
    img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    
    # 4. Generate HTML report
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Functional Distance Heatmap</title>
        <style>
            body {{ font-family: sans-serif; margin: 20px; }}
            h1 {{ color: #333; }}
            .heatmap-container {{ border: 1px solid #ccc; padding: 20px; display: inline-block; }}
        </style>
    </head>
    <body>
        <h1>Pairwise Functional Distance Heatmap</h1>
        <p>This heatmap visualizes the functional distance between all gene variants, computed using the 'functional_distance' metric from the Evo2 scoring report. Lower distances (darker colors) indicate higher functional similarity.</p>
        
        <div class="heatmap-container">
            <img src="data:image/png;base64,{img_base64}" alt="Functional Distance Heatmap" style="max-width: 100%; height: auto;">
        </div>
    </body>
    </html>
    """
    
    # Save the HTML report
    with open(HEATMAP_HTML_OUTPUT_PATH, 'w') as f:
        f.write(html_content)
        
    print(f"Successfully generated and saved HTML heatmap report to: {HEATMAP_HTML_OUTPUT_PATH}")

except ImportError:
    print("\nWarning: Libraries (matplotlib/seaborn) required for visualization were not found. Skipping heatmap generation.")
except Exception as e:
    print(f"\nCritical error during heatmap generation or HTML saving: {e}")

print("\n--- Task B Complete ---")
