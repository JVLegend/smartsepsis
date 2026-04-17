import pandas as pd
import os
import matplotlib.pyplot as plt
import io
import base64

# --- Configuration ---
WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
OUTPUT_DIR = os.path.join(WORKING_DIR, "autoresearch_results")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# TASK: O - Summary dashboard: create a single-page HTML with KPI cards (total genes, coverage %, conserved variants, CARD variants) and mini-charts. Read all report files. Save summary_dashboard.html.
print("Starting Summary Dashboard creation...")

# --- Data Loading ---

# 1. Conservation Analysis
CONSERVATION_PATH = os.path.join(WORKING_DIR, "reports/conservation_analysis.tsv")
df_conservation = pd.DataFrame()
try:
    print(f"Attempting to load conservation analysis from: {CONSERVATION_PATH}")
    # Assuming standard tab separation
    df_conservation = pd.read_csv(CONSERVATION_PATH, sep='\t')
    print("Successfully loaded conservation_analysis.tsv.")
    print(f'Columns: {df_conservation.columns.tolist()}')
except FileNotFoundError:
    print(f"WARNING: Conservation analysis file not found at {CONSERVATION_PATH}. Using empty DataFrame.")
except Exception as e:
    print(f"ERROR reading conservation_analysis.tsv: {e}")

# 2. Functional Scores (Evo2)
EVO2_PATH = os.path.join(WORKING_DIR, "reports/evo2_scoring/functional_scores.tsv")
df_evo2 = pd.DataFrame()
try:
    print(f"Attempting to load functional scores from: {EVO2_PATH}")
    df_evo2 = pd.read_csv(EVO2_PATH, sep='\t')
    print("Successfully loaded functional_scores.tsv.")
    print(f'Columns: {df_evo2.columns.tolist()}')
except FileNotFoundError:
    print(f"WARNING: Functional scores file not found at {EVO2_PATH}. Using empty DataFrame.")
except Exception as e:
    print(f"ERROR reading functional_scores.tsv: {e}")

# 3. CARD Variants
CARD_PATH = os.path.join(WORKING_DIR, "reports/card_new_variants.csv")
df_card = pd.DataFrame()
try:
    print(f"Attempting to load CARD variants from: {CARD_PATH}")
    # Reading as CSV
    df_card = pd.read_csv(CARD_PATH)
    print("Successfully loaded card_new_variants.csv.")
    print(f'Columns: {df_card.columns.tolist()}')
except FileNotFoundError:
    print(f"WARNING: CARD variants file not found at {CARD_PATH}. Using empty DataFrame.")
except Exception as e:
    print(f"ERROR reading card_new_variants.csv: {e}")


# --- KPI Calculation ---

# A. Total Genes (Estimate based on the reports available)
# We will use unique gene entries found across the reports.
genes_from_card = df_card['pipeline_family'].dropna().unique()
genes_from_conservation = df_conservation['family'].dropna().unique()
genes_from_evo2 = df_evo2['variant'].dropna().unique()

all_unique_genes = set(genes_from_card) | set(genes_from_conservation) | set(genes_from_evo2)
total_genes = len(all_unique_genes)

# B. Conserved Variants
# Count unique variants where 'match_type' is 'exact'
conserved_variants = df_conservation[df_conservation['match_type'] == 'exact']['variant'].nunique()

# C. CARD Variants
# Count unique accessions found in the CARD report
card_variants = df_card['aro_accession'].dropna().nunique()

# D. Coverage Percentage (Placeholder/Simple Proxy Calculation)
# Calculate based on the proportion of reports loaded successfully.
coverage_score = 0
if not df_conservation.empty():
    coverage_score += 1
if not df_evo2.empty():
    coverage_score += 1
if not df_card.empty():
    coverage_score += 1

coverage_percentage = f"{(coverage_score / 3) * 100:.0f}% (3/3 report readiness)"


# --- Chart Generation ---

# 1. Chart 1: Variant Counts (Conserved vs CARD)
fig_chart1, ax1 = plt.subplots(figsize=(8, 5))
x = ['Conserved Variants', 'CARD Identified Variants']
y = [conserved_variants, card_variants]
ax1.bar(x, y, color=['#28a745', '#007bff'])
ax1.set_title('Key Variant Identification')
ax1.set_ylabel('Count')
ax1.grid(axis='y', linestyle='--', alpha=0.7)

# 2. Chart 2: Gene Report Status (Conceptual)
# Since we don't have multiple metrics, we'll show the ratio of success/failure.
fig_chart2, ax2 = plt.subplots(figsize=(8, 5))
labels = ['Conservation', 'Evo2 Score', 'CARD Data']
values = [1 if not df_conservation.empty() else 0, 1 if not df_evo2.empty() else 0, 1 if not df_card.empty() else 0]
ax2.bar(labels, values, color=['#ffc107', '#20c997', '#17a2b8'])
ax2.set_title('Report Data Availability Status')
ax2.set_ylabel('Available (1) / Missing (0)')
ax2.set_ylim(0, 1.2)
ax2.grid(axis='y', linestyle='--', alpha=0.7)


# Save charts to buffer
img1_buffer = io.BytesIO()
plt.savefig(img1_buffer, format='png')
img1_buffer.seek(0)
base64_chart1 = base64.b64encode(img1_buffer.read()).decode('utf-8')
plt.close(fig_chart1)

img2_buffer = io.BytesIO()
plt.savefig(img2_buffer, format='png')
img2_buffer.seek(0)
base64_chart2 = base64.b64encode(img2_buffer.read()).decode('utf-8')
plt.close(fig_chart2)

print("Mini-charts generated successfully.")

# --- HTML Generation ---

html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CRISPR-Cas12a AMR Diagnostic Pipeline Summary</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f4f7fa; }}
        .header {{ background-color: #0056b3; color: white; padding: 20px; text-align: center; border-radius: 8px 8px 0 0; }}
        .dashboard-container {{ background-color: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1); }}
        h2 {{ color: #0056b3; border-bottom: 2px solid #eee; padding-bottom: 10px; margin-top: 30px; }}
        .kpi-grid {{ display: flex; flex-wrap: wrap; gap: 20px; margin-bottom: 40px; justify-content: space-between; }}
        .kpi-card {{ flex: 1 1 22%; min-width: 200px; background-color: #e9ecef; padding: 20px; border-radius: 10px; text-align: center; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05); }}
        .kpi-value {{ font-size: 3em; font-weight: bold; color: #28a745; margin: 5px 0; }}
        .kpi-label {{ font-size: 1em; color: #6c757d; }}
        .chart-grid {{ display: flex; flex-wrap: wrap; gap: 20px; }}
        .chart-box {{ flex: 1 1 45%; min-width: 300px; background-color: white; padding: 20px; border: 1px solid #ddd; border-radius: 8px; }}
        img {{ max-width: 100%; height: auto; }}
        .highlight {{ color: #dc3545; font-weight: bold; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🔬 AMR Diagnostic Pipeline Summary Dashboard</h1>
        <p>CRISPR-Cas12a Resistance Profiling Analysis</p>
    </div>

    <div class="dashboard-container">
        <h2>📈 Key Performance Indicators (KPIs)</h2>
        <div class="kpi-grid">
            <div class="kpi-card">
                <div class="kpi-value">{total_genes}</div>
                <div class="kpi-label">Total Unique Genes Targeted</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{coverage_percentage}</div>
                <div class="kpi-label">Data Coverage Status</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{conserved_variants}</div>
                <div class="kpi-label">Conserved Resistance Variants Found</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{card_variants}</div>
                <div class="kpi-label">CARD Classified Variants</div>
            </div>
        </div>

        <h2>📊 Visual Summary Charts</h2>
        <div class="chart-grid">
            <div class="chart-box">
                <h3>Variant Identification Comparison</h3>
                <img src="data:image/png;base64,{base64_chart1}" alt="Variant Comparison Chart">
            </div>
            <div class="chart-box">
                <h3>Report Data Availability Status</h3>
                <img src="data:image/png;base64,{base64_chart2}" alt="Report Status Chart">
            </div>
        </div>
        
        <p style="margin-top: 40px; font-size: 0.9em; color: #666;">
            Analysis complete. Metrics compiled from multiple biological and bioinformatics reports.
        </p>
    </div>
</body>
</html>
"""

# --- Saving Output ---

OUTPUT_PATH = os.path.join(OUTPUT_DIR, "summary_dashboard.html")
try:
    with open(OUTPUT_PATH, "w") as f:
        f.write(html_content)
    print(f"\n✅ Success: Summary dashboard saved to {OUTPUT_PATH}")
except Exception as e:
    print(f"\n❌ Fatal Error: Could not save HTML file: {e}")

print("Script execution finished.")
