# TASK: Variant coverage heatmap: read conservation_analysis.tsv, create a gene-family vs match-type (exact/near/none) heatmap. Use matplotlib base64 in HTML. Save coverage_heatmap.html.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import base64
import io
import os

def create_coverage_heatmap():
    """
    Loads conservation data, generates a heatmap showing the count of variants 
    (coverage) based on Gene Family vs. Match Type, and saves it as a self-contained HTML file.
    """
    
    # --- Setup Paths ---
    WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
    INPUT_FILE = os.path.join(WORKING_DIR, "reports/conservation_analysis.tsv")
    OUTPUT_DIR = os.path.join(WORKING_DIR, "autoresearch_results/")
    OUTPUT_FILE = os.path.join(OUTPUT_DIR, "coverage_heatmap.html")

    print("--- Starting Coverage Heatmap Generation ---")

    # Ensure output directory exists
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        print(f"Output directory ensured: {OUTPUT_DIR}")
    except OSError as e:
        print(f"Error creating output directory {OUTPUT_DIR}: {e}")
        return

    # 1. Load Data
    try:
        print(f"1. Loading data from: {INPUT_FILE}")
        df = pd.read_csv(INPUT_FILE, sep='\t')
        print(f"   Data loaded successfully. Shape: {df.shape}")
    except FileNotFoundError:
        print(f"ERROR: Input file not found at {INPUT_FILE}. Exiting.")
        return
    except Exception as e:
        print(f"ERROR loading data: {e}. Exiting.")
        return

    # 2. Data Preprocessing and Aggregation
    # We count the occurrences of (family, match_type) pairs
    print("2. Aggregating variant counts by Gene Family and Match Type...")
    
    # Use the size() method on a groupby object to count rows
    heatmap_data_series = df.groupby(['family', 'match']).size()
    
    # Create a DataFrame suitable for heatmap visualization
    heatmap_df = heatmap_data_series.unstack(fill_value=0)

    # 3. Generate Heatmap Visualization
    print("3. Generating heatmap visualization...")
    
    plt.figure(figsize=(10, 8))
    
    # Use seaborn to create the heatmap
    # Annotating with the count value is crucial for interpretation
    sns.heatmap(
        heatmap_df, 
        annot=True, 
        fmt='d', # format as integer
        cmap="viridis", # Nice color map
        linewidths=.5, 
        linecolor='lightgray',
        cbar_kws={'label': 'Number of Variants/Coverage'}
    )
    
    plt.title('Variant Coverage Heatmap: Gene Family vs. Match Type', fontsize=16)
    plt.ylabel('Gene Family', fontsize=12)
    plt.xlabel('Match Type (Exact/Near/None)', fontsize=12)
    plt.tight_layout()

    # 4. Save Plot to BytesIO Buffer
    print("4. Saving plot to in-memory buffer (PNG)...")
    
    # Save the figure into a BytesIO object
    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format='png', bbox_inches='tight')
    plt.close() # Close the plot to free memory

    # Get the base64 encoding of the image bytes
    image_base64 = base64.b64encode(img_buffer.getvalue()).decode('utf-8')
    img_src = f"data:image/png;base64,{image_base64}"

    # 5. Generate and Save HTML Report
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Variant Coverage Heatmap</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1 {{ color: #333; border-bottom: 2px solid #eee; padding-bottom: 10px; }}
            .container {{ max-width: 1200px; margin: auto; }}
            img {{ max-width: 100%; height: auto; border: 1px solid #ddd; padding: 10px; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>CRISPR-Cas12a Variant Coverage Heatmap</h1>
            <p>This heatmap visualizes the count of observed variants (coverage) across different <b>Gene Families</b> (y-axis) segmented by their <b>Match Type</b> (x-axis). The intensity of the color represents the frequency of variants found.</p>
            
            <h2>Data Visualization</h2>
            <img src="{img_src}" alt="Variant Coverage Heatmap">

            <h3>Interpretation</h3>
            <ul>
                <li><strong>Color Intensity:</strong> Indicates the frequency (count) of variants found for that specific Gene Family and Match Type combination.</li>
                <li><strong>Gene Family:</strong> The grouping of related bacterial resistance genes.</li>
                <li><strong>Match Type:</strong> Indicates the level of conservation/similarity: 
                    <ul>
                        <li><strong>Exact:</strong> Perfect sequence match observed.</li>
                        <li><strong>Near:</strong> Sequence similarity with a few mismatches.</li>
                        <li><strong>None:</strong> No detectable sequence homology (lowest coverage).</li>
                    </ul>
                </li>
            </ul>
        </div>
    </body>
    </html>
    """
    
    try:
        with open(OUTPUT_FILE, "w") as f:
            f.write(html_content)
        print(f"\n✅ Success! Coverage heatmap saved to: {OUTPUT_FILE}")
    except Exception as e:
        print(f"\n❌ ERROR saving HTML file: {e}")

if __name__ == "__main__":
    create_coverage_heatmap()
