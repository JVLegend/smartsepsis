import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import base64
import io
import os
from typing import List

def generate_coverage_heatmap(base_dir: str, output_dir: str):
    """
    Reads the conservation_analysis.tsv, creates a heatmap of Gene Family vs 
    Match Type coverage, and saves the result to an HTML file.
    """
    input_file = os.path.join(base_dir, 'reports/conservation_analysis.tsv')
    output_html = os.path.join(output_dir, 'coverage_heatmap.html')

    print("--- Task H: Variant Coverage Heatmap Generation ---")

    # 1. Read the conservation_analysis.tsv
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"Successfully loaded conservation_analysis.tsv from: {input_file}")
        print(f"Columns: {df.columns.tolist()}")
    except FileNotFoundError:
        print(f"ERROR: Input file not found at {input_file}. Exiting.")
        return
    except Exception as e:
        print(f"ERROR reading file: {e}. Exiting.")
        return

    # 2. Data Preparation: Count occurrences of (family, match_type)
    # We need to count how many records exist for each combination.
    
    # Ensure 'family' and 'match_type' columns exist before proceeding
    required_cols = ['family', 'match_type']
    if not all(col in df.columns for col in required_cols):
        print(f"ERROR: The DataFrame must contain {required_cols}. Found: {df.columns.tolist()}")
        return

    coverage_df = pd.crosstab(df['family'], df['match_type'])
    print("\nGenerated coverage matrix (Rows: Family, Columns: Match Type).")

    # 3. Create and customize the heatmap
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Use the count values for the heatmap
    im = ax.imshow(coverage_df.values, cmap='viridis', aspect='auto')
    
    # Add colorbar
    plt.colorbar(im, ax=ax, label='Observed Variant Count')
    
    # Set labels
    ax.set_xticks(np.arange(len(coverage_df.columns)))
    ax.set_yticks(np.arange(len(coverage_df.index)))
    
    # Set column labels (Match Type)
    plt.xlabel("Match Type (Column)")
    plt.title("Gene Family Coverage Heatmap (Variant Count by Match Type)")
    
    # Set row labels (Family)
    ax.set_xticklabels(coverage_df.columns, rotation=45, ha="right")
    ax.set_yticklabels(coverage_df.index)

    # Annotate cells with count numbers
    for i in range(len(coverage_df.index)):
        for j in range(len(coverage_df.columns)):
            count = coverage_df.iloc[i, j]
            plt.text(j, i, str(count), ha='center', va='center', color='black', fontsize=9)

    # Adjust layout for better visualization
    plt.tight_layout()
    
    # 4. Save plot to base64
    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format='png', bbox_inches='tight')
    plt.close(fig)
    
    base64_img = base64.b64encode(img_buffer.getvalue()).decode('utf-8')
    
    # 5. Create HTML report
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Variant Coverage Heatmap</title>
        <style>
            body {{ font-family: sans-serif; margin: 20px; }}
            h1 {{ color: #333; }}
            img {{ max-width: 100%; height: auto; border: 1px solid #ccc; padding: 10px; }}
        </style>
    </head>
    <body>
        <h1>Gene Family vs. Match Type Coverage Heatmap</h1>
        <p>This heatmap visualizes the count of observed variants (N) categorized by the associated Gene Family (Y-axis) and the match type of the sequence (X-axis).</p>
        <p>Darker/Brighter colors indicate higher variant coverage.</p>
        <img src="data:image/png;base64,{base64_img}" alt="Coverage Heatmap">
        <p><em>Source Data: conservation_analysis.tsv</em></p>
    </body>
    </html>
    """

    # 6. Save the HTML file
    try:
        with open(output_html, 'w') as f:
            f.write(html_content)
        print(f"\n[SUCCESS] Coverage Heatmap saved successfully to: {output_html}")
    except Exception as e:
        print(f"[ERROR] Failed to write output HTML file: {e}")


if __name__ == "__main__":
    # Define absolute paths as required
    BASE_WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
    OUTPUT_RESULTS_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/autoresearch_results/"

    # Ensure the output directory exists
    os.makedirs(OUTPUT_RESULTS_DIR, exist_ok=True)

    generate_coverage_heatmap(BASE_WORKING_DIR, OUTPUT_RESULTS_DIR)
