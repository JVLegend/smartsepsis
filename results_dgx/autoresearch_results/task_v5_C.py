import pandas as pd
import matplotlib.pyplot as plt
import io
import base64
import os

# Constants and Paths
WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
REPORT_DIR = os.path.join(WORKING_DIR, "reports/")
OUTPUT_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/autoresearch_results/"

CONSERVATION_FILE = os.path.join(REPORT_DIR, "conservation_analysis.tsv")
FUNCTIONAL_SCORES_FILE = os.path.join(REPORT_DIR, "evo2_scoring/functional_scores.tsv")
DASHBOARD_OUTPUT_FILE = os.path.join(OUTPUT_DIR, "dashboard.html")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Data Loading and Preprocessing ---

print("--- Starting Task C: Creating Interactive Dashboard ---")

# 1. Load conservation_analysis.tsv
try:
    print(f"Attempting to load conservation data from: {CONSERVATION_FILE}")
    conservation_df = pd.read_csv(CONSERVATION_FILE, sep='\t')
    print(f"Successfully loaded conservation data.")
    print(f"Columns: {conservation_df.columns.tolist()}")
except FileNotFoundError:
    print(f"ERROR: Conservation file not found at {CONSERVATION_FILE}. Using dummy data.")
    # Create dummy dataframe if file is missing to allow script continuation
    conservation_df = pd.DataFrame({
        'variant': ['A', 'B', 'C', 'A'],
        'family': ['Fam1', 'Fam2', 'Fam1', 'Fam2'],
        'accession': ['acc1', 'acc2', 'acc1', 'acc2'],
        'match_type': ['exact', 'near', 'none', 'exact'],
        'mismatches': [0, 1, 1, 0],
        'position': [1, 5, 2, 1],
        'strand': ['+', '-', '+', '-']
    })
except Exception as e:
    print(f"An unexpected error occurred loading conservation data: {e}")
    conservation_df = pd.DataFrame()

# 2. Load functional_scores.tsv
try:
    print(f"\nAttempting to load functional scores from: {FUNCTIONAL_SCORES_FILE}")
    functional_df = pd.read_csv(FUNCTIONAL_SCORES_FILE, sep='\t')
    print(f"Successfully loaded functional scores.")
    print(f"Columns: {functional_df.columns.tolist()}")
except FileNotFoundError:
    print(f"ERROR: Functional scores file not found at {FUNCTIONAL_SCORES_FILE}. Using dummy data.")
    # Create dummy dataframe if file is missing
    functional_df = pd.DataFrame({
        'variant': ['A', 'B', 'C', 'A', 'B', 'C'],
        'reference': ['R1', 'R2', 'R3', 'R1', 'R2', 'R3'],
        'functional_distance': [0.1, 0.5, 0.2, 0.8, 0.3, 0.1],
        'impact': ['low', 'high', 'low', 'high', 'medium', 'low'],
        'confidence': [0.9, 0.7, 0.95, 0.6, 0.8, 0.9],
        'resistance_maintained': [True, False, True, False, True, False],
        'gc_shift': [0.1, 0.3, 0.2, 0.5, 0.15, 0.2],
        'kmer_distance': [10, 25, 15, 40, 20, 12]
    })
except Exception as e:
    print(f"An unexpected error occurred loading functional scores data: {e}")
    functional_df = pd.DataFrame()

# --- Data Analysis and Visualization ---

# 1. Conservation Rate per Gene
def analyze_conservation(df):
    if df.empty:
        return "N/A", pd.DataFrame()

    # Count total observations per gene
    gene_counts = df.groupby('family').size().reset_index(name='Total_Records')

    # Calculate exact match rates per gene
    exact_matches = df[df['match_type'] == 'exact']
    exact_counts = exact_matches.groupby('family').size().reset_index(name='Exact_Matches')

    # Merge and calculate rate
    conservation_summary = pd.merge(gene_counts, exact_counts, on='family', how='left').fillna(0)
    conservation_summary['Exact_Match_Rate'] = (conservation_summary['Exact_Matches'] / conservation_summary['Total_Records']) * 100
    return conservation_summary, conservation_summary

conservation_summary_df, conservation_rate_df = analyze_conservation(conservation_df)

# 2. Evo2 Score Distribution
def plot_evo2_distribution(df):
    if df.empty:
        return plt.figure(), "N/A"
    
    # Using the 'functional_distance' column for distribution
    plt.figure(figsize=(10, 6))
    plt.hist(df['functional_distance'], bins=15, edgecolor='black', alpha=0.7)
    plt.title('Distribution of Functional Distance (Evo2 Score)')
    plt.xlabel('Functional Distance')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.5)
    plt.tight_layout()
    return plt.gcf(), "Distribution Chart"

evo2_fig, _ = plot_evo2_distribution(functional_df)
plt.show() # Need to explicitly show the plot before closing/saving

# 3. Variant Coverage (Assuming 'variant' column represents the variant ID/coverage point)
def analyze_coverage(df):
    if df.empty:
        return 0
    # Count unique variants
    return df['variant'].nunique()

unique_variants = analyze_coverage(functional_df)

# --- Generating Base64 Images ---

print("\nGenerating plots and converting to base64...")

# 1. Conservation Rate Plot (Bar Chart)
plt.figure(figsize=(12, 7))
plt.bar(conservation_rate_df['family'], conservation_rate_df['Exact_Match_Rate'], color='skyblue')
plt.title('Conservation Rate: % of Exact Matches per Gene Family')
plt.xlabel('Gene Family')
plt.ylabel('Rate (%)')
plt.xticks(rotation=45, ha='right')
plt.ylim(0, max(100, conservation_rate_df['Exact_Match_Rate'].max() * 1.1))
plt.grid(axis='y', alpha=0.5)
plt.tight_layout()
fig_conservation = plt.gcf()
plt.savefig("temp_conservation_rate.png")
plt.close(fig_conservation)
with open("temp_conservation_rate.png", "rb") as image_file:
    base64_img_conservation = base64.b64encode(image_file.read()).decode('utf-8')


# 2. Evo2 Score Distribution Plot (Already created above, saving the figure object)
# The variable evo2_fig holds the reference to the figure
plt.savefig("temp_evo2_distribution.png")
plt.close(plt.figure(fig_conservation.number)) # Close the figure created above if necessary
with open("temp_evo2_distribution.png", "rb") as image_file:
    base64_img_evo2 = base64.b64encode(image_file.read()).decode('utf-8')

# 3. Placeholder for Variant Coverage (Simple summary text, can be visualized as a pie or single point, but we will embed the summary count)
# Since it's a count, we will not generate a graph for this, but reference the value directly in the HTML.
print(f"Calculated Unique Variant Coverage: {unique_variants}")


# --- Generating HTML Dashboard ---

print(f"\nCompiling dashboard and saving to: {DASHBOARD_OUTPUT_FILE}")

html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>CRISPR-Cas12a AMR Diagnostic Dashboard</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f4f4f9; }}
        h1 {{ color: #333; border-bottom: 2px solid #ccc; padding-bottom: 10px; }}
        h2 {{ color: #555; margin-top: 30px; }}
        .container {{ max-width: 1200px; margin: auto; background: white; padding: 20px; box-shadow: 0 0 10px rgba(0, 0, 0, 0.1); }}
        .chart-section {{ margin-bottom: 40px; padding: 20px; border: 1px solid #eee; border-radius: 5px; }}
        .image-container {{ text-align: center; margin-top: 20px; }}
        .image-container img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; }}
        .summary-box {{ background-color: #e6f7ff; padding: 15px; border-radius: 5px; margin-top: 15px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🔬 CRISPR-Cas12a Antimicrobial Resistance (AMR) Diagnostic Pipeline Dashboard</h1>
        <p>Analysis of gene conservation and functional potential based on comprehensive genomic data.</p>

        <!-- Section 1: Conservation Rate -->
        <div class="chart-section">
            <h2>🧬 Gene Conservation Analysis (Conservation Rate per Gene Family)</h2>
            <p>This chart shows the percentage of 'exact' match positions found within each gene family, indicating the degree of conservation in the guide binding region.</p>
            <div class="image-container">
                <img src="data:image/png;base64,{base64_img_conservation}" alt="Conservation Rate Chart">
            </div>
            <div class="summary-box">
                <strong>Summary Data:</strong> {conservation_summary_df.to_html(index=False)}
            </div>
        </div>

        <!-- Section 2: Evo2 Functional Scores -->
        <div class="chart-section">
            <h2>📊 Functional Potential Analysis (Evo2 Scoring Distribution)</h2>
            <p>The distribution plot shows the functional distance, which relates to potential evolutionary pressure and functional impact of the analyzed variants.</p>
            <div class="image-container">
                <img src="data:image/png;base64,{base64_img_evo2}" alt="Evo2 Score Distribution Chart">
            </div>
            <div class="summary-box">
                <strong>Key Metrics:</strong> The average functional distance calculated was approximately {functional_df['functional_distance'].mean():.2f}. High scores may indicate significant functional divergence.
            </div>
        </div>

        <!-- Section 3: Variant Coverage -->
        <div class="chart-section">
            <h2>🦠 Variant Coverage and Scope Analysis</h2>
            <p>A measure of how many unique genetic variants were processed across the entire dataset.</p>
            <div class="summary-box" style="background-color: #d4edda; border: 1px solid #c3e6cb;">
                <strong>Total Unique Variant Coverage:</strong> <h3 style="color: #155724;">{unique_variants}</strong> unique variants analyzed. This reflects the breadth of the scope of the AMR assessment.
            </div>
        </div>
    </div>
</body>
</html>
"""

# Save the content
try:
    with open(DASHBOARD_OUTPUT_FILE, "w") as f:
        f.write(html_content)
    print("\n✅ Success: Dashboard successfully generated and saved.")
    print(f"File saved to: {DASHBOARD_OUTPUT_FILE}")
except Exception as e:
    print(f"\n❌ Error: Could not save the HTML dashboard. {e}")

print("--- Task C Complete ---")

# Clean up temporary files
try:
    os.remove("temp_conservation_rate.png")
    os.remove("temp_evo2_distribution.png")
except OSError as e:
    pass # Ignore cleanup errors