import pandas as pd
import matplotlib.pyplot as plt
import io
import base64
import os
from datetime import datetime

# --- Configuration ---
BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
REPORT_DIR = os.path.join(BASE_DIR, "reports/")
EVO2_SCORES_PATH = os.path.join(REPORT_DIR, "evo2_scoring/functional_scores.tsv")
CONSERVATION_PATH = os.path.join(REPORT_DIR, "conservation_analysis.tsv")
OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results/")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "dashboard.html")

# TASK: CRISPR-Cas12a AMR diagnostic pipeline dashboard generation

def load_data(filepath: str, sheet_name: str) -> pd.DataFrame | None:
    """Loads data with error handling."""
    print(f"-> Attempting to load data from: {filepath}")
    try:
        # Assuming all specified reports are tab-separated (TSV)
        df = pd.read_csv(filepath, sep='\t')
        print(f"   Successfully loaded data. Shape: {df.shape}")
        return df
    except FileNotFoundError:
        print(f"   ERROR: File not found at {filepath}. Cannot generate dashboard.")
        return None
    except Exception as e:
        print(f"   ERROR loading {filepath}: {e}")
        return None

def generate_conservation_chart(df_cons: pd.DataFrame) -> tuple[str, str]:
    """
    Generates a bar chart showing average match counts per gene family.
    Returns (base64_image_data, title).
    """
    print("-> Generating Conservation Rate Plot...")
    
    # Calculate average match count and total variants per family
    conservation_summary = df_cons.groupby('family')[['match']].mean().reset_index()
    conservation_summary.columns = ['Family', 'Average Matches']
    
    plt.figure(figsize=(12, 6))
    plt.bar(conservation_summary['Family'], conservation_summary['Average Matches'], color='teal')
    plt.xlabel("Gene Family")
    plt.ylabel("Average Number of Matches")
    plt.title("Average Match Rate per Gene Family")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    # Save plot to bytes buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    plt.close()
    
    base64_image = base64.b64encode(buf.getvalue()).decode('utf-8')
    return base64_image, "Gene Family Conservation Rates"

def generate_evo2_score_chart(df_evo2: pd.DataFrame) -> tuple[str, str]:
    """
    Generates a histogram distribution chart for Functional Distance (Impact).
    Returns (base64_image_data, title).
    """
    print("-> Generating Evo2 Score Distribution Plot...")
    
    plt.figure(figsize=(10, 6))
    
    # Use a histogram to show the distribution of functional distance scores
    plt.hist(df_evo2['functional_distance'].dropna(), bins=30, color='darkorange', edgecolor='black')
    
    plt.xlabel("Functional Distance Score (Impact)")
    plt.ylabel("Frequency of Variants")
    plt.title("Distribution of Functional Distance Scores (Evo2)")
    plt.grid(axis='y', alpha=0.7)
    plt.tight_layout()
    
    # Save plot to bytes buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    plt.close()
    
    base64_image = base64.b64encode(buf.getvalue()).decode('utf-8')
    return base64_image, "Evo2 Functional Distance Distribution"

def create_dashboard_html(
    base64_cons: str, title_cons: str,
    base64_evo2: str, title_evo2: str,
    data_cons: pd.DataFrame | None,
    data_evo2: pd.DataFrame | None
) -> str:
    """Constructs the full HTML dashboard content."""
    
    # Calculate summary metrics for display
    summary_cons = ""
    if data_cons is not None:
        avg_matches = data_cons['match'].mean()
        unique_families = data_cons['family'].nunique()
        summary_cons = f"**Summary:** Found data for {unique_families} unique gene families. The average match count across all analyzed sites is **{avg_matches:.2f}**."
    
    summary_evo2 = ""
    if data_evo2 is not None:
        mean_distance = data_evo2['functional_distance'].mean()
        std_distance = data_evo2['functional_distance'].std()
        summary_evo2 = f"**Summary:** The average functional distance score is **{mean_distance:.2f}**, with a standard deviation of **{std_distance:.2f}**. This indicates moderate variability in predicted functional impact."


    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>CRISPR-Cas12a AMR Diagnostic Dashboard</title>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">
    <style>
        body {{ font-family: Arial, sans-serif; padding: 20px; background-color: #f8f9fa; }}
        .card {{ margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1 {{ color: #0056b3; border-bottom: 2px solid #0056b3; padding-bottom: 10px; }}
        h2 {{ color: #333; margin-top: 20px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1 class="text-center">CRISPR-Cas12a AMR Diagnostic Pipeline Results</h1>
        <p class="lead text-center">Analysis Dashboard created on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

        <!-- SECTION 1: Conservation Analysis -->
        <div class="card">
            <div class="card-header bg-primary text-white">
                <h2>📈 {title_cons}</h2>
            </div>
            <div class="card-body">
                <p>{summary_cons}</p>
                <div class="row text-center mb-4">
                    <div class="col-md-6">
                        <img src="data:image/png;base64,{base64_cons}" class="img-fluid border rounded" alt="{title_cons} chart">
                    </div>
                </div>
                <h3>Raw Data Snapshot (Top 5 Families):</h3>
                <pre class="bg-light p-3 border overflow-auto" style="max-height: 300px;">{data_cons.head().to_markdown(index=False)}</pre>
            </div>
        </div>

        <!-- SECTION 2: Functional Scoring Analysis -->
        <div class="card">
            <div class="card-header bg-success text-white">
                <h2>📊 {title_evo2}</h2>
            </div>
            <div class="card-body">
                <p>{summary_evo2}</p>
                <div class="row text-center mb-4">
                    <div class="col-md-12">
                        <img src="data:image/png;base64,{base64_evo2}" class="img-fluid border rounded" alt="{title_evo2} chart">
                    </div>
                </div>
                <h3>Raw Data Snapshot (Top 5 Scores):</h3>
                <pre class="bg-light p-3 border overflow-auto" style="max-height: 300px;">{data_evo2.head().to_markdown(index=False)}</pre>
            </div>
        </div>

        <footer class="text-center mt-5 pt-3 border-top">
            <p class="text-muted">Dashboard generated by expert Python pipeline. All analysis based on provided reports.</p>
        </footer>
    </div>
</body>
</html>
"""
    return html_content

def main():
    """Main execution function to run the pipeline."""
    print("==================================================")
    print("STARTING AMR DIAGNOSTIC DASHBOARD PIPELINE")
    print("==================================================")

    # 1. Setup Directory
    try:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        print(f"[*] Output directory ensured: {OUTPUT_DIR}")
    except Exception as e:
        print(f"CRITICAL ERROR: Could not create output directory {OUTPUT_DIR}. {e}")
        return

    # 2. Load Data
    df_cons = load_data(CONSERVATION_PATH, "conservation_analysis.tsv")
    df_evo2 = load_data(EVO2_SCORES_PATH, "functional_scores.tsv")

    if df_cons is None or df_evo2 is None:
        print("\n--- FAILURE ---")
        print("Cannot proceed with dashboard generation due to missing or corrupted input files.")
        return

    # 3. Generate Visualizations and Base64 Images
    
    # Conservation Plot
    base64_cons, title_cons = generate_conservation_chart(df_cons)
    
    # Evo2 Plot
    base64_evo2, title_evo2 = generate_evo2_score_chart(df_evo2)

    # 4. Generate HTML Content
    print("\n-> Compiling HTML dashboard...")
    html_output = create_dashboard_html(
        base64_cons, title_cons,
        base64_evo2, title_evo2,
        df_cons, df_evo2
    )

    # 5. Save Output
    try:
        with open(OUTPUT_FILE, 'w') as f:
            f.write(html_output)
        print("\n==================================================")
        print(f"✅ SUCCESS! Dashboard saved to: {OUTPUT_FILE}")
        print("==================================================")
    except Exception as e:
        print(f"CRITICAL ERROR: Failed to write dashboard to {OUTPUT_FILE}. {e}")

if __name__ == "__main__":
    # Ensure matplotlib can handle required backend operations
    # This is useful if running in environments like Jupyter notebooks where plt.close() behavior is critical.
    plt.rcParams['figure.dpi'] = 150
    main()
