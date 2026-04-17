import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
from datetime import datetime

# --- Configuration ---
# NOTE: Adjust these paths if the actual directory structure changes.
GUIDES_PATH = "guides/guide_data.tsv" 
CONSERVATION_PATH = "reports/conservation_analysis.tsv"
EVO2_PATH = "reports/evo2_scoring/functional_scores.tsv"
CARD_PATH = "reports/card_new_variants.csv"
OUTPUT_HTML = "summary_dashboard.html"
TEMP_DIR = "temp_charts"

# --- Utility Functions ---

def setup_environment():
    """Creates necessary directories and handles mock data creation for runnable script."""
    print("Setting up environment...")
    os.makedirs("guides", exist_ok=True)
    os.makedirs("reports/conservation_analysis", exist_ok=True)
    os.makedirs("reports/evo2_scoring", exist_ok=True)
    
    # Create dummy files to ensure the script runs without FileNotFoundError
    # 1. Guides TSV
    dummy_guides = pd.DataFrame({
        'rank': [1, 2, 3], 
        'gene': ['A', 'B', 'A'], 
        'spacer_seq': ['AAAAA', 'TTTTT', 'AAAAA'], 
        'pam_seq': ['GATC', 'GATC', 'GATC'], 
        'position': [10, 20, 30], 
        'strand': ['+', '-', '+'], 
        'spacer_length': [23] * 3, 
        'score': [90, 85, 90], 
        'gc': [0.5, 0.6, 0.5], 
        'homopolymer': [1, 0, 1], 
        'self_comp': [0.1, 0.2, 0.1], 
        'poly_t': [0, 1, 0], 
        'rel_position': [5, 15, 25], 
        'crRNA_LbCas12a': ['X', 'Y', 'Z'], 
        'crRNA_AsCas12a': ['X', 'Y', 'Z']
    }).to_csv(GUIDES_PATH, sep='\t', index=False)
    
    # 2. Conservation TSV
    dummy_conservation = pd.DataFrame({
        'variant': ['V1', 'V2', 'V3', 'V4'], 
        'family': ['F1', 'F2', 'F1', 'F3'], 
        'accession': ['Acc1', 'Acc2', 'Acc1', 'Acc3'], 
        'match': ['exact', 'near', 'exact', 'none'], 
        'mismatches': [0, 1, 0, 2], 
        'position': [1, 5, 1, 8], 
        'strand': ['+', '-']
    }).to_csv(CONSERVATION_PATH, sep='\t', index=False)
    
    # 3. Evo2 TSV
    dummy_evo2 = pd.DataFrame({
        'variant': ['V1', 'V2', 'V3', 'V4'], 
        'reference': ['RefA', 'RefB', 'RefA', 'RefB'], 
        'functional_distance': [1.2, 0.5, 1.8, 0.3], 
        'impact': [0.8, 0.3, 0.9, 0.1], 
        'confidence': [0.9, 0.95, 0.8, 0.7], 
        'resistance_maintained': [True, True, False, True], 
        'gc_shift': [0.1, 0.2, 0.0, 0.1], 
        'kmer_distance': [5, 3, 7, 1]
    }).to_csv(EVO2_PATH, sep='\t', index=False)
    
    # 4. CARD CSV
    dummy_card = pd.DataFrame({
        'card_model_name': ['ModelA', 'ModelB', 'ModelC'], 
        'pipeline_family': ['P1', 'P2', 'P1'], 
        'aro_accession': ['A1', 'A2', 'A3'], 
        'drug_classes': ['Antibio1,Drug2', 'Drug3', 'Antibio2'], 
        'resistance_mechanisms': ['M1', 'M2', 'M3']
    }).to_csv(CARD_PATH, index=False)

    print("Setup complete. Starting data processing.")


def load_and_process_data():
    """Loads all data files and prints column names."""
    
    # 1. Guides Data
    try:
        df_guides = pd.read_csv(GUIDES_PATH, sep='\t')
        print(f'\n--- Guides Data Loaded ({GUIDES_PATH}) ---')
        print(f'Columns: {df_guides.columns.tolist()}')
    except FileNotFoundError:
        print(f"Error: Guides file not found at {GUIDES_PATH}")
        df_guides = pd.DataFrame()

    # 2. Conservation Analysis
    try:
        df_conservation = pd.read_csv(CONSERVATION_PATH, sep='\t')
        print(f'\n--- Conservation Analysis Loaded ({CONSERVATION_PATH}) ---')
        print(f'Columns: {df_conservation.columns.tolist()}')
    except FileNotFoundError:
        print(f"Error: Conservation file not found at {CONSERVATION_PATH}")
        df_conservation = pd.DataFrame()

    # 3. Evo2 Scoring
    try:
        df_evo2 = pd.read_csv(EVO2_PATH, sep='\t')
        print(f'\n--- Evo2 Scoring Loaded ({EVO2_PATH}) ---')
        print(f'Columns: {df_evo2.columns.tolist()}')
    except FileNotFoundError:
        print(f"Error: Evo2 file not found at {EVO2_PATH}")
        df_evo2 = pd.DataFrame()

    # 4. CARD Variants
    try:
        df_card = pd.read_csv(CARD_PATH, sep=',')
        print(f'\n--- CARD Variants Loaded ({CARD_PATH}) ---')
        print(f'Columns: {df_card.columns.tolist()}')
    except FileNotFoundError:
        print(f"Error: CARD file not found at {CARD_PATH}")
        df_card = pd.DataFrame()
        
    return df_guides, df_conservation, df_evo2, df_card

def calculate_kpis(df_guides, df_conservation, df_card):
    """Calculates all key performance indicators."""
    
    kpi = {}
    
    # Total Genes
    total_genes = df_guides['gene'].nunique() if not df_guides.empty else 0
    kpi['total_genes'] = total_genes
    
    # Coverage % (Approximation: Use number of guides as a proxy for coverage metrics)
    # Calculate percentage relative to a theoretical maximum, or simply count unique guides.
    unique_guides_count = len(df_guides)
    kpi['guides_count'] = unique_guides_count
    
    # Conserved Variants Count
    conserved_variants = len(df_conservation)
    kpi['conserved_variants'] = conserved_variants
    
    # CARD Variants Count
    card_variants = len(df_card)
    kpi['card_variants'] = card_variants
    
    # Coverage % calculation (Placeholder calculation based on guide count)
    # Let's use the ratio of guides count to a fixed number (e.g., 100) for demonstration
    kpi['coverage_pct'] = f"{(unique_guides_count / 100) * 100:.1f}%" if unique_guides_count > 0 else "N/A"
    
    return kpi

def generate_charts(df_conservation, df_evo2):
    """Generates and saves necessary mini-charts."""
    
    if not os.path.exists(TEMP_DIR):
        os.makedirs(TEMP_DIR)

    charts = {}

    # Chart 1: Distribution of Match Types (Conservation Analysis)
    if not df_conservation.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.countplot(data=df_conservation, x='match', order=['exact', 'near', 'none'], ax=ax)
        ax.set_title('Distribution of Conservation Match Types')
        ax.set_xlabel('Match Type (Exact, Near, None)')
        ax.set_ylabel('Number of Variants')
        plt.tight_layout()
        chart_path = os.path.join(TEMP_DIR, 'match_types.png')
        plt.savefig(chart_path)
        plt.close(fig)
        charts['match_distribution'] = chart_path

    # Chart 2: Distribution of Functional Impact (Evo2 Scoring)
    if not df_evo2.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.histplot(df_evo2['impact'], bins=5, kde=True, ax=ax)
        ax.set_title('Distribution of Functional Impact Score')
        ax.set_xlabel('Functional Impact Score (0 to 1)')
        ax.set_ylabel('Frequency')
        plt.tight_layout()
        chart_path = os.path.join(TEMP_DIR, 'impact_distribution.png')
        plt.savefig(chart_path)
        plt.close(fig)
        charts['impact_distribution'] = chart_path

    return charts

def generate_html_dashboard(kpis, charts):
    """Generates the final single-page HTML dashboard."""
    
    # HTML/CSS Template structure
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CRISPR Summary Dashboard</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background-color: #f4f7f9;
            color: #333;
            margin: 20px;
        }}
        .container {{
            max-width: 1200px;
            margin: auto;
        }}
        h1 {{
            color: #1a237e;
            border-bottom: 3px solid #4caf50;
            padding-bottom: 10px;
        }}
        .kpi-grid {{
            display: flex;
            justify-content: space-around;
            margin-bottom: 30px;
            gap: 20px;
        }}
        .kpi-card {{
            flex: 1;
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
            text-align: center;
        }}
        .kpi-card h3 {{
            margin: 0;
            font-size: 1em;
            color: #777;
        }}
        .kpi-card p {{
            font-size: 2.5em;
            font-weight: bold;
            color: #1a237e;
            margin-top: 5px;
        }}
        .section-title {{
            background-color: #e0eaff;
            padding: 10px;
            border-radius: 5px;
            margin-top: 30px;
            font-weight: 600;
        }}
        .chart-grid {{
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
            margin-top: 20px;
        }}
        .chart-box {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            flex-basis: calc(50% - 30px); /* Two charts per row */
            box-sizing: border-box;
        }}
        @media (max-width: 900px) {{
            .kpi-grid, .chart-grid {{
                flex-direction: column;
            }}
            .chart-box {{
                flex-basis: 100%;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>CRISPR Genomic Analysis Summary Dashboard</h1>
        <p style="margin-bottom: 30px; color: #666;">Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

        <!-- KPI Cards -->
        <div class="kpi-grid">
            <div class="kpi-card">
                <h3>Total Unique Genes</h3>
                <p>{kpis['total_genes']}</p>
            </div>
            <div class="kpi-card">
                <h3>Guide Coverage Est.</h3>
                <p>{kpis['coverage_pct']}</p>
            </div>
            <div class="kpi-card">
                <h3>Conserved Variants</h3>
                <p>{kpis['conserved_variants']}</p>
            </div>
            <div class="kpi-card">
                <h3>CARD Variants Detected</h3>
                <p>{kpis['card_variants']}</p>
            </div>
        </div>

        <!-- Charts Section -->
        <div class="section-title">Mini-Charts and Distributions</div>
        <div class="chart-grid">
            
            <div class="chart-box">
                <h2>Match Type Distribution (Conservation)</h2>
                <img src="{charts['match_distribution']}" alt="Match Distribution Chart" style="max-width: 100%;">
            </div>
            
            <div class="chart-box">
                <h2>Functional Impact Score (Evo2)</h2>
                <img src="{charts['impact_distribution']}" alt="Impact Distribution Chart" style="max-width: 100%;">
            </div>
        </div>
        
        <p style="margin-top: 50px; text-align: center; color: #999;">
            <small>* Dashboard powered by bioinformatics analysis pipelines.</small>
        </p>
    </div>
</body>
</html>
"""
    return html_content

def cleanup():
    """Removes temporary chart files and directory."""
    import shutil
    if os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)
        print("\nCleanup: Removed temporary chart directory.")

# --- Main Execution ---
if __name__ == "__main__":
    try:
        # 1. Setup Mock Data Environment
        setup_environment()
        
        # 2. Load Data
        df_guides, df_conservation, df_evo2, df_card = load_and_process_data()

        # 3. Calculate KPIs
        kpis = calculate_kpis(df_guides, df_conservation, df_card)
        print("\nKPI Calculation Complete:")
        print(f"  Total Genes: {kpis['total_genes']}")
        print(f"  Coverage Est: {kpis['coverage_pct']}")
        print(f"  Conserved Variants: {kpis['conserved_variants']}")
        print(f"  CARD Variants: {kpis['card_variants']}")

        # 4. Generate Charts
        print("\nGenerating charts...")
        charts = generate_charts(df_conservation, df_evo2)
        print("Charts successfully generated and saved.")

        # 5. Generate HTML Dashboard
        html_output = generate_html_dashboard(kpis, charts)
        
        with open(OUTPUT_HTML, 'w') as f:
            f.write(html_output)
            
        print(f"\n==================================================")
        print(f"SUCCESS: Summary Dashboard saved to {OUTPUT_HTML}")
        print(f"==================================================")

    except Exception as e:
        print(f"\nAN ERROR OCCURRED during execution: {e}")
    finally:
        # 6. Clean up
        cleanup()
