import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import base64
import io
import os
import sys

# Set up global paths
WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
OUTPUT_DIR = os.path.join(WORKING_DIR, "autoresearch_results")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Helper Functions ---

def create_base64_image(fig: plt.Figure) -> str:
    """Converts a matplotlib figure object to a Base64 encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches='tight')
    plt.close(fig)
    img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    return f"data:image/png;base64,{img_base64}"

def load_data_with_robustness(file_path: str, file_type: str) -> pd.DataFrame | None:
    """Loads a file using try/except and handles various file extensions."""
    print(f"Loading {file_type} data from: {file_path}")
    try:
        if file_type == 'tsv':
            df = pd.read_csv(file_path, sep='\t')
        elif file_type == 'csv':
            df = pd.read_csv(file_path)
        else:
            print(f"Warning: Unknown file type for {file_path}. Skipping.")
            return None
        return df
    except FileNotFoundError:
        print(f"!!! ERROR: File not found at {file_path}. Skipping this section.")
        return None
    except pd.errors.EmptyDataError:
        print(f"!!! ERROR: The file {file_path} is empty. Skipping this section.")
        return None
    except Exception as e:
        print(f"!!! CRITICAL ERROR loading {file_path}: {e}. Skipping this section.")
        return None

def generate_dashboard_report() -> str:
    """
    Reads, processes, visualizes, and generates the complete HTML report.
    """
    print("\n--- Starting Data Loading ---")

    # 1. Load Reports Data
    # reports/conservation_analysis.tsv
    conservation_path = os.path.join(WORKING_DIR, "reports/conservation_analysis.tsv")
    df_conservation = load_data_with_robustness(conservation_path, 'tsv')

    # reports/evo2_scoring/functional_scores.tsv
    evo2_path = os.path.join(WORKING_DIR, "reports/evo2_scoring/functional_scores.tsv")
    df_evo2 = load_data_with_robustness(evo2_path, 'tsv')

    # reports/card_new_variants.csv
    card_path = os.path.join(WORKING_DIR, "reports/card_new_variants.csv")
    df_card = load_data_with_robustness(card_path, 'csv')

    # reports/covariance_probes/*.tsv (Assume we process all in this directory)
    covariance_dir = os.path.join(WORKING_DIR, "reports/covariance_probes/")
    covariance_files = [os.path.join(covariance_dir, f) for f in os.listdir(covariance_dir) if f.endswith('.tsv')]
    all_covariance_dfs = []
    print("\n--- Loading Covariance Probes ---")
    for i, f_path in enumerate(covariance_files):
        df_cov = load_data_with_robustness(f_path, 'tsv')
        if df_cov is not None:
            all_covariance_dfs.append(df_cov)
    
    if all_covariance_dfs:
        # Concatenate all covariance probes for a single summary
        df_covariance = pd.concat(all_covariance_dfs, ignore_index=True)
    else:
        df_covariance = pd.DataFrame()

    print("\n--- Starting Visualization and Report Generation ---")

    # --- Section 1: Conservation Analysis ---
    conservation_html = ""
    if df_conservation is not None and not df_conservation.empty:
        try:
            # Simple Example Visualization: Distribution of Mismatches
            fig_con = sns.countplot(y='family', x='mismatches', data=df_conservation, order=df_conservation['family'].value_counts().index)
            con_base64 = create_base64_image(fig_con.figure_)
            
            conservation_html = f"""
            <h3>🧬 Conservation Analysis</h3>
            <p>Visualization of observed mutation diversity (Mismatches) across different gene families.</p>
            <div style="text-align: center; margin: 20px 0;">
                <img src="{con_base64}" alt="Mismatches by Family" style="max-width: 100%; height: auto;">
            </div>
            <h4>Key Metrics Summary</h4>
            <div style="max-width: 800px; margin: 20px 0; overflow-x: auto;">
                {df_conservation.head(10).to_html(classes='table table-striped', index=False)}
            </div>
            """
        except Exception as e:
            conservation_html = f"<h3>🧬 Conservation Analysis</h3><p style='color: red;'>Error generating conservation plot or data: {e}</p>"


    # --- Section 2: Evo2 Scoring ---
    evo2_html = ""
    if df_evo2 is not None and not df_evo2.empty:
        try:
            # Example Visualization: Correlation between Functional Distance and Impact
            fig_evo2 = sns.scatterplot(x='functional_distance', y='impact', data=df_evo2)
            evo2_base64 = create_base64_image(fig_evo2.figure_)
            
            evo2_html = f"""
            <h3>🧬 Evo2 Functional Scoring</h3>
            <p>Assessment of evolutionary distance and predicted functional impact of variants.</p>
            <div style="text-align: center; margin: 20px 0;">
                <img src="{evo2_base64}" alt="Evo2 Distance vs Impact" style="max-width: 100%; height: auto;">
            </div>
            <h4>Summary of Predicted Changes</h4>
            <div style="max-width: 800px; margin: 20px 0; overflow-x: auto;">
                {df_evo2.head(10).to_html(classes='table table-striped', index=False)}
            </div>
            """
        except Exception as e:
            evo2_html = f"<h3>🧬 Evo2 Functional Scoring</h3><p style='color: red;'>Error generating Evo2 plot or data: {e}</p>"


    # --- Section 3: CARD Model ---
    card_html = ""
    if df_card is not None and not df_card.empty:
        try:
            # Example Visualization: Unique resistance mechanisms count (simple bar chart simulation)
            mechanisms = df_card['resistance_mechanisms'].unique()
            if len(mechanisms) > 0:
                # Create dummy data for visualization if columns are complex
                data_for_plot = pd.Series(mechanisms)
                fig_card = sns.countplot(y=data_for_plot.str.split(', ').str.extract(r'(\w+)'), order=list(data_for_plot.str.split(', ').str.extract(r'(\w+)').value_counts().index))
                # Since multi-label plotting is complex, we will just plot unique count vs index
                plt.figure(figsize=(10, 4))
                plt.barh([f"Mechanism {i+1}" for i in range(len(mechanisms))], [1] * len(mechanisms))
                plt.yticks(range(len(mechanisms)), mechanisms)
                plt.xlabel("Count")
                plt.title("Unique Resistance Mechanisms Detected")
                card_base64 = create_base64_image(plt.gcf())
                
                card_html = f"""
                <h3>🛡️ CARD Model Prediction</h3>
                <p>Identification of potential drug resistance classes and mechanisms based on sequence homology.</p>
                <div style="text-align: center; margin: 20px 0;">
                    <img src="{card_base64}" alt="Resistance Mechanisms" style="max-width: 100%; height: auto;">
                </div>
                <h4>Top Resistance Profiles Found</h4>
                <div style="max-width: 800px; margin: 20px 0; overflow-x: auto;">
                    {df_card[['card_model_name', 'drug_classes', 'resistance_mechanisms']].head(10).to_html(classes='table table-striped', index=False)}
                </div>
                """
            else:
                card_html = "<h3>🛡️ CARD Model Prediction</h3><p>No resistance mechanism data found.</p>"
        except Exception as e:
            card_html = f"<h3>🛡️ CARD Model Prediction</h3><p style='color: red;'>Error generating CARD plot or data: {e}</p>"


    # --- Section 4: Biophysical Covariance Probes ---
    covariance_html = ""
    if not df_covariance.empty:
        try:
            # Example Visualization: Average GC shift vs. Kmer Distance
            # Using a simple summary plot since the data structure is assumed to be various probes
            fig_cov = sns.scatterplot(x='kmer_distance', y='gc_shift', data=df_covariance.dropna())
            cov_base64 = create_base64_image(fig_cov.figure_)
            
            covariance_html = f"""
            <h3>🔬 Biophysical Covariance Probes</h3>
            <p>Summary of structural and physicochemical changes predicted by various probes.</p>
            <div style="text-align: center; margin: 20px 0;">
                <img src="{cov_base64}" alt="GC Shift vs Kmer Distance" style="max-width: 100%; height: auto;">
            </div>
            <h4>Overall Feature Summary</h4>
            <div style="max-width: 800px; margin: 20px 0; overflow-x: auto;">
                {df_covariance[['gc_shift', 'kmer_distance', 'impact']].describe().to_html(classes='table table-striped', index=True)}
            </div>
            """
        except Exception as e:
            covariance_html = f"<h3>🔬 Biophysical Covariance Probes</h3><p style='color: red;'>Error generating covariance plot or data: {e}</p>"
    else:
        covariance_html = """
        <h3>🔬 Biophysical Covariance Probes</h3>
        <p>No covariance probe data was available for analysis.</p>
        """


    # --- Final HTML Structure Assembly ---
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>CRISPR-Cas12a AMR Diagnostic Pipeline Report</title>
        <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">
        <style>
            body {{ font-family: Arial, sans-serif; padding: 20px; background-color: #f4f7f9; }}
            .container-fluid {{ max-width: 1200px; }}
            .section-header {{ padding: 20px; background-color: #e9ecef; border-radius: 8px; margin-bottom: 30px; }}
            h1 {{ color: #0056b3; margin-top: 0; }}
            h3 {{ color: #343a40; border-bottom: 2px solid #dee2e6; padding-bottom: 5px; margin-top: 20px; }}
        </style>
    </head>
    <body>
        <div class="container-fluid">
            <header class="text-center mb-5">
                <h1 style="color: #28a745;">🧬 CRISPR-Cas12a AMR Diagnostic Pipeline Report</h1>
                <p class="lead">Comprehensive Analysis of Antimicrobial Resistance (AMR) Variants</p>
                <p class="text-muted">Generated from conservation, structural, and functional prediction models.</p>
            </header>

            <!-- Section 1: Conservation -->
            <section id="conservation" class="mb-5">
                <div class="section-header">
                    <h2>🦠 Core Genetic Conservation Analysis</h2>
                </div>
                {conservation_html}
            </section>

            <!-- Section 2: Evo2 Scoring -->
            <section id="evo2" class="mb-5">
                <div class="section-header">
                    <h2>📈 Evolutionary and Functional Scoring (Evo2)</h2>
                </div>
                {evo2_html}
            </section>

            <!-- Section 3: CARD Model -->
            <section id="card" class="mb-5">
                <div class="section-header">
                    <h2>🛡️ Drug Resistance Model Prediction (CARD)</h2>
                </div>
                {card_html}
            </section>

            <!-- Section 4: Covariance Probes -->
            <section id="covariance" class="mb-5">
                <div class="section-header">
                    <h2>🔬 Biophysical Covariance Probes</h2>
                </div>
                {covariance_html}
            </section>

        </div>
    </body>
    </html>
    """
    return html_content

# --- Main Execution ---

if __name__ == "__main__":
    try:
        print("Starting the comprehensive report generation process...")
        html_report_content = generate_dashboard_report()
        
        output_file_path = os.path.join(WORKING_DIR, "autoresearch_results", "comprehensive_report.html")
        
        with open(output_file_path, "w") as f:
            f.write(html_report_content)
        
        print("\n===================================================")
        print(f"✅ Success! Comprehensive HTML report saved to:")
        print(f"   {output_file_path}")
        print("===================================================")

    except Exception as e:
        print(f"\n🛑 A critical error occurred during the report generation: {e}", file=sys.stderr)
        sys.exit(1)
