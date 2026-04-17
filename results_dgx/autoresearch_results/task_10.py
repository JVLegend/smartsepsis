import pandas as pd
import matplotlib.pyplot as plt
import base64
import io
import os
import sys

def generate_card_enrichment_report():
    """
    Reads the CARD variant reports, performs enrichment counting (pipeline_family vs drug_class),
    and generates a stacked bar chart saved as an HTML file with embedded base64 plot.
    """
    
    # --- Configuration ---
    WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
    INPUT_FILE = os.path.join(WORKING_DIR, "reports/card_new_variants.csv")
    OUTPUT_DIR = os.path.join(WORKING_DIR, "autoresearch_results/")
    OUTPUT_HTML_FILE = os.path.join(OUTPUT_DIR, "card_enrichment.html")

    print(f"TASK: J - CARD variant enrichment analysis starting.")
    print(f"1. Reading input file: {INPUT_FILE}")

    # 1. Data Loading and Preparation
    try:
        df = pd.read_csv(INPUT_FILE)
        print("   [SUCCESS] Data loaded successfully.")
    except FileNotFoundError:
        print(f"   [ERROR] Input file not found at {INPUT_FILE}. Exiting.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print("   [ERROR] Input file is empty. Exiting.")
        sys.exit(1)
    except Exception as e:
        print(f"   [ERROR] An unexpected error occurred during file reading: {e}. Exiting.")
        sys.exit(1)

    # Check for required columns
    required_columns = ['pipeline_family', 'drug_classes']
    if not all(col in df.columns for col in required_columns):
        print(f"   [ERROR] Missing required columns. Expected: {required_columns}. Found: {df.columns.tolist()}. Exiting.")
        sys.exit(1)

    # 2. Data Transformation: Exploding drug_classes for accurate counting
    print("2. Processing data: Splitting and counting variants by drug class...")
    
    try:
        # Split the comma-separated drug classes and explode the list
        df['drug_class'] = df['drug_classes'].str.split(',').explode()
        
        # Group and count the occurrences
        counts = df.groupby(['pipeline_family', 'drug_class']).size().reset_index(name='Count')
        print(f"   [SUCCESS] Data aggregated. Unique pipeline families: {counts['pipeline_family'].nunique()}.")

    except Exception as e:
        print(f"   [ERROR] Error during data transformation/grouping: {e}. Exiting.")
        sys.exit(1)

    # 3. Visualization Preparation (Pivoting)
    print("3. Preparing data structure for stacked bar chart...")
    
    # Use pivot to get a count for every drug class within every pipeline family
    pivot_df = counts.pivot(index='pipeline_family', columns='drug_class', values='Count').fillna(0)
    
    # Ensure all drug classes are represented, even if they start at zero
    drug_classes = sorted(pivot_df.columns)
    
    # 4. Plot Generation
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 8))
    
    bottom = None
    
    # Create stacked bars iteratively
    for i, drug_class in enumerate(drug_classes):
        counts_series = pivot_df[drug_class].values
        plt.bar(
            pivot_df.index, 
            counts_series, 
            label=str(drug_class), 
            bottom=bottom
        )
        if bottom is None:
            bottom = counts_series
        else:
            bottom = bottom + counts_series

    ax.set_xlabel("Pipeline Family", fontsize=14)
    ax.set_ylabel("Number of Variants", fontsize=14)
    ax.set_title("Enrichment of Variants by Drug Class Across Pipeline Families", fontsize=16, pad=20)
    ax.legend(title="Drug Class", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout for legend
    
    # 5. Saving the Figure (Base64 Encoding)
    print("4. Encoding plot and saving report...")

    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Save the plot to a BytesIO object in PNG format
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    plt.close(fig)

    # Encode the bytes object to Base64
    base64_encoded_image = base64.b64encode(buf.getvalue()).decode('utf-8')
    
    # 6. Generating the HTML File
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CARD Variant Enrichment Analysis</title>
    <style>
        body {{ font-family: Arial, sans-serif; padding: 20px; background-color: #f4f4f9; }}
        .container {{ background-color: #ffffff; padding: 30px; border-radius: 10px; box-shadow: 0 4px 8px rgba(0,0,0,0.1); }}
        h1 {{ color: #333; border-bottom: 2px solid #007bff; padding-bottom: 10px; }}
        .plot-container {{ text-align: center; margin-top: 20px; }}
        img {{ max-width: 100%; height: auto; display: block; }}
        p.metadata {{ margin-top: 20px; font-style: italic; color: #666; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>CRISPR-Cas12a AMR Diagnostic Pipeline: Drug Class Enrichment</h1>
        <p>This stacked bar chart visualizes the count of reported antimicrobial resistance (AMR) variants, grouped by the initial screening pipeline family, showing how these variants are distributed across different drug classes.</p>
        
        <div class="plot-container">
            <img src="data:image/png;base64,{base64_encoded_image}" alt="Stacked Bar Chart of AMR Variant Enrichment">
        </div>
        
        <p class="metadata">Data Source: reports/card_new_variants.csv</p>
    </div>
</body>
</html>
"""

    try:
        with open(OUTPUT_HTML_FILE, 'w') as f:
            f.write(html_content)
        print(f"   [SUCCESS] Report saved to {OUTPUT_HTML_FILE}")
    except Exception as e:
        print(f"   [CRITICAL ERROR] Could not save the HTML file: {e}")

    print("TASK: J - Analysis pipeline finished.")


if __name__ == "__main__":
    generate_card_enrichment_report()
