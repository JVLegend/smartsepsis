import os
import pandas as pd
import matplotlib.pyplot as plt
import base64
import io
import numpy as np
from datetime import datetime

# --- Configuration ---
BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results/")
REPORT_DIR = os.path.join(BASE_DIR, "reports/")

# Input Paths
GUIDES_PATH = os.path.join(BASE_DIR, "guides/*.tsv")
CONSERVATION_PATH = os.path.join(BASE_DIR, "reports/conservation_analysis.tsv")
EVO2_SCORE_PATH = os.path.join(BASE_DIR, "reports/evo2_scoring/functional_scores.tsv")
CARD_PATH = os.path.join(BASE_DIR, "reports/card_new_variants.csv")
COV_PROBES_PATH = os.path.join(BASE_DIR, "reports/covariance_probes/")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Target output directory set to: {OUTPUT_DIR}")

# --- Helper Functions ---

def load_tsv(file_pattern, required_cols):
    """Loads a TSV file, handling multiple files and missing columns."""
    try:
        all_files = [f for f in os.listdir(os.path.join(BASE_DIR, '')) if f.startswith("guides") and f.endswith(".tsv")]
        if not all_files:
            print("Warning: No FASTA/TSV files found in guides/ folder.")
            return pd.DataFrame()

        df_list = []
        print(f"Loading multiple guide files (up to {len(all_files)})...")
        for filename in all_files:
            filepath = os.path.join(BASE_DIR, filename)
            try:
                df = pd.read_csv(filepath, sep='\t')
                df_list.append(df)
            except Exception as e:
                print(f"Error loading {filepath}: {e}")
        
        if not df_list:
             print("Error: Could not load any guide data.")
             return pd.DataFrame()

        combined_df = pd.concat(df_list, ignore_index=True)
        # Ensure required columns are present
        for col in required_cols:
            if col not in combined_df.columns:
                print(f"Warning: Required column '{col}' missing from combined dataset.")
        return combined_df
    except Exception as e:
        print(f"CRITICAL ERROR loading TSV data: {e}")
        return pd.DataFrame()

def load_csv(file_path):
    """Loads a single CSV file."""
    try:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found at {file_path}")
        return pd.read_csv(file_path)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        return pd.DataFrame()
    except Exception as e:
        print(f"CRITICAL ERROR loading CSV data from {file_path}: {e}")
        return pd.DataFrame()

def generate_base64_chart(fig):
    """Saves a matplotlib figure to BytesIO and encodes it to base64."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches='tight')
    plt.close(fig)
    return base64.b64encode(buf.getbuffer()).decode("ascii")

# --- Data Processing & KPI Calculation ---

print("\n--- Starting Data Processing ---")

# 1. Load Required Datasets
guides_df = load_tsv(GUIDES_PATH, ['guide_id', 'sequence', 'gene_family', 'pam_site'])
conservation_df = load_csv(CONSERVATION_PATH)
evo2_df = load_csv(EVO2_SCORE_PATH)
card_df = load_csv(CARD_PATH)

# 2. KPI Calculations
total_genes = guides_df['guide_id'].nunique() if not guides_df.empty else 0

# Coverage % (Simulating calculation based on the presence of functional scores)
# Assuming a maximum reference coverage (e.g., 1000) and calculating average detected functional distance.
avg_functional_distance = evo2_df['functional_distance'].mean() if not evo2_df.empty and 'functional_distance' in evo2_df.columns else 0.0
# Simplified KPI: Using the average functional distance as a proxy for coverage/detection depth
coverage_percent = min(100.0, avg_functional_distance * 10) # Example scaling

# Conserved Variants Count
conserved_variants_count = conservation_df.shape[0] if not conservation_df.empty else 0

# CARD Variants Count
card_variants_count = card_df.shape[0] if not card_df.empty else 0

# 3. Chart Data Generation

# Chart 1: Resistance Mechanism Distribution (Pie Chart)
resistance_mechanisms = card_df['resistance_mechanisms'].dropna().str.split(',').explode()
resistance_counts = resistance_mechanisms.value_counts().head(10)
plt.figure(figsize=(8, 8))
plt.pie(resistance_counts.values, labels=resistance_counts.index[:5].tolist(), autopct='%1.1f%%', startangle=140, explode=[0.05] * len(resistance_counts))
plt.title('Top 5 Dominant Resistance Mechanisms')
fig_chart1 = plt.gcf()
base64_chart1 = generate_base64_chart(fig_chart1)
print("Generated Chart 1: Resistance Mechanisms.")


# Chart 2: Variant Impact Distribution (Bar Chart)
impact_counts = evo2_df['impact'].value_counts().to_frame('Count')
plt.figure(figsize=(10, 6))
impact_counts.plot(kind='bar', color='skyblue')
plt.title('Distribution of Functional Impact')
plt.xlabel('Impact Category')
plt.ylabel('Number of Variants')
plt.xticks(rotation=45, ha='right')
plt.grid(axis='y', alpha=0.5)
fig_chart2 = plt.gcf()
base64_chart2 = generate_base64_chart(fig_chart2)
print("Generated Chart 2: Impact Distribution.")


# --- HTML Dashboard Generation ---

# Format metrics
kpis = {
    "Total Genes Identified": f"{total_genes:,}",
    "Avg. Functional Coverage %": f"{coverage_percent:.2f}%",
    "Conserved Variants Found": f"{conserved_variants_count:,}",
    "CARD Known Variants": f"{card_variants_count:,}"
}

# HTML Structure
html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CRISPR-Cas12a AMR Diagnostic Dashboard</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f4f7fa; }}
        .container {{ max-width: 1200px; margin: auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1); }}
        h1 {{ color: #0056b3; border-bottom: 2px solid #0056b3; padding-bottom: 10px; }}
        .kpi-grid {{ display: flex; justify-content: space-around; flex-wrap: wrap; margin-bottom: 30px; }}
        .kpi-card {{ background: #e6f0ff; border-left: 5px solid #007bff; padding: 15px 20px; margin: 10px; border-radius: 6px; flex-basis: calc(25% - 30px); min-width: 200px; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05); }}
        .kpi-card h3 {{ margin: 0; font-size: 1.1em; color: #555; }}
        .kpi-card p {{ font-size: 2.2em; font-weight: bold; color: #007bff; margin-top: 5px; }}
        .chart-section {{ margin-top: 30px; padding: 20px; border: 1px solid #ddd; border-radius: 8px; }}
        h2 {{ color: #333; border-bottom: 1px solid #eee; padding-bottom: 10px; }}
        .chart-grid {{ display: flex; justify-content: space-around; flex-wrap: wrap; gap: 20px; }}
        .chart-card {{ flex: 1 1 45%; min-width: 350px; text-align: center; }}
        .chart-card img {{ max-width: 100%; height: auto; border: 1px solid #ccc; border-radius: 4px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 CRISPR-Cas12a AMR Diagnostic Summary Dashboard</h1>
        <p>Generated on: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>

        <div class="kpi-grid">
            <div class="kpi-card">
                <h3>Total Genes Identified</h3>
                <p>{kpis['Total Genes Identified']}</p>
            </div>
            <div class="kpi-card">
                <h3>Avg. Functional Coverage</h3>
                <p>{kpis['Avg. Functional Coverage %']}</p>
            </div>
            <div class="kpi-card">
                <h3>Conserved Variants Found</h3>
                <p>{kpis['Conserved Variants Found']}</p>
            </div>
            <div class="kpi-card">
                <h3>CARD Known Variants</h3>
                <p>{kpis['CARD Known Variants']}</p>
            </div>
        </div>

        <div class="chart-section">
            <h2>🔬 Key Distribution Charts</h2>
            <div class="chart-grid">
                <div class="chart-card">
                    <h3>Resistance Mechanism Distribution</h3>
                    <img src="data:image/png;base64,{base64_chart1}" alt="Resistance Mechanism Chart">
                </div>
                <div class="chart-card">
                    <h3>Functional Impact Distribution</h3>
                    <img src="data:image/png;base64,{base64_chart2}" alt="Impact Distribution Chart">
                </div>
            </div>
        </div>
    </div>
</body>
</html>
"""

# Inject variables into the template
final_html = html_template.format(
    kpis=kpis,
    base64_chart1=base64_chart1,
    base64_chart2=base64_chart2
)

# 4. Save Output
output_filepath = os.path.join(OUTPUT_DIR, "summary_dashboard.html")

try:
    with open(output_filepath, "w") as f:
        f.write(final_html)
    print("\n" + "="*50)
    print("SUCCESS: Summary dashboard created.")
    print(f"File saved to: {output_filepath}")
    print("="*50)
except Exception as e:
    print(f"\nERROR: Could not save the final HTML dashboard. {e}")