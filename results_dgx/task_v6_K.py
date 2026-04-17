import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from io import BytesIO
import base64

# --- Configuration and Setup ---
BASE_DIR = ""
GUIDES_DIR = "guides/"
REPORT_DIR = "reports/"
EVO2_DIR = "reports/evo2_scoring/"

# --- Helper Functions ---

def load_tsv(filepath, sep='\t'):
    """Loads a TSV file using pandas."""
    try:
        df = pd.read_csv(filepath, sep=sep)
        print(f"Successfully loaded {filepath}")
        print(f'Columns: {df.columns.tolist()}')
        return df
    except FileNotFoundError:
        print(f"Warning: File not found at {filepath}")
        return pd.DataFrame()
    except pd.errors.EmptyDataError:
        print(f"Warning: Empty file at {filepath}")
        return pd.DataFrame()

def load_csv(filepath, sep=','):
    """Loads a CSV file using pandas."""
    try:
        df = pd.read_csv(filepath, sep=sep)
        print(f"Successfully loaded {filepath}")
        print(f'Columns: {df.columns.tolist()}')
        return df
    except FileNotFoundError:
        print(f"Warning: File not found at {filepath}")
        return pd.DataFrame()
    except pd.errors.EmptyDataError:
        print(f"Warning: Empty file at {filepath}")
        return pd.DataFrame()

def create_base64_chart(fig):
    """Saves a matplotlib figure to a base64 PNG string."""
    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight')
    plt.close(fig)
    image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    return f'<img src="data:image/png;base64,{image_base64}" alt="Chart">'

# --- 1. Data Loading ---

print("--- Starting Data Loading ---")

# A. Load Guides Data (Iterating through files in guides/)
guides_data = []
guide_files = [f for f in os.listdir(GUIDES_DIR) if f.endswith('.tsv')]

print(f"\n--- Loading {len(guide_files)} Guide Files ---")
for filename in guide_files:
    filepath = os.path.join(GUIDES_DIR, filename)
    df_guide = load_tsv(filepath)
    if not df_guide.empty:
        guides_data.append(df_guide)

df_guides = pd.concat(guides_data, ignore_index=True)
print("\nGuides Data Consolidation Complete.")


# B. Load Conservation Data
CONSERVATION_FILE = os.path.join(REPORT_DIR, 'conservation_analysis.tsv')
df_conservation = load_tsv(CONSERVATION_FILE)

# C. Load Evo2 Data
EVO2_FILE = os.path.join(EVO2_DIR, 'functional_scores.tsv')
df_evo2 = load_tsv(EVO2_FILE)

# D. Load CARD Data
CARD_FILE = os.path.join(REPORT_DIR, 'card_new_variants.csv')
df_card = load_csv(CARD_FILE)

# --- 2. Data Processing & Visualization ---

print("\n--- Generating Charts and Reports ---")

# 2.1 Conservation Chart (Match Counts)
con_counts = df_conservation['match'].value_counts()
plt.figure(figsize=(8, 5))
sns.barplot(x=con_counts.index, y=con_counts.values)
plt.title('Distribution of Conservation Matches')
plt.xlabel('Match Type (exact, near, none)')
plt.ylabel('Count')
plt.xticks(rotation=0)
chart_conservation = create_base64_chart(plt.figure(figsize=(8, 5)))


# 2.2 Evo2 Chart (Impact Distribution)
impact_counts = df_evo2['impact'].value_counts()
plt.figure(figsize=(8, 5))
sns.countplot(x=impact_counts.index, data=df_evo2)
plt.title('Distribution of Functional Impact Scores')
plt.xlabel('Impact')
plt.ylabel('Frequency')
chart_evo2 = create_base64_chart(plt.figure(figsize=(8, 5)))


# --- 3. HTML Generation ---

def create_html_dashboard(df_guides, df_conservation, df_evo2, df_card, chart_conservation, chart_evo2):
    """Assembles all data and charts into a comprehensive HTML structure."""

    # HTML for the Guides Section
    guides_html = f"""
    <h2>🧬 CRISPR Guide Data Summary</h2>
    <p>This section summarizes all loaded CRISPR guide data, providing details on spacer sequences, gene association, and genomic context.</p>
    <p>Total Records: {len(df_guides)}</p>
    <div style="max-width: 100%; overflow-x: auto;">
        {df_guides.to_html(classes='table table-striped', index=False, border=0)}
    </div>
    """

    # HTML for the Conservation Section
    conservation_html = f"""
    <h2>🛡️ Conservation Analysis Report</h2>
    <p>Analysis of sequence conservation based on reported variants, identifying patterns of match types (exact, near, none).</p>
    <h3>Match Type Distribution</h3>
    <div style="text-align: center; margin: 20px 0;">
        {chart_conservation}
    </div>
    <h3>Raw Data Table</h3>
    <div style="max-width: 100%; overflow-x: auto;">
        {df_conservation.to_html(classes='table table-striped', index=False, border=0)}
    </div>
    """

    # HTML for the Evo2 Section
    evo2_html = f"""
    <h2>⚡ Evo2 Functional Scoring</h2>
    <p>Assessment of functional scores (e.g., impact, distance) for variants, indicating potential changes in protein function or structural stability.</p>
    <h3>Impact Score Distribution</h3>
    <div style="text-align: center; margin: 20px 0;">
        {chart_evo2}
    </div>
    <h3>Raw Data Table</h3>
    <div style="max-width: 100%; overflow-x: auto;">
        {df_evo2.to_html(classes='table table-striped', index=False, border=0)}
    </div>
    """

    # HTML for the CARD Section
    card_html = f"""
    <h2>🦠 CARD Antibiotic Resistance Data</h2>
    <p>Identification of potential antibiotic resistance mechanisms associated with new variants, drawn from the CARD database.</p>
    <h3>Summary of Variants and Resistance Profiles</h3>
    <div style="max-width: 100%; overflow-x: auto;">
        {df_card.to_html(classes='table table-striped', index=False, border=0)}
    </div>
    """

    # Combine all sections into the final structure
    html_dashboard = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Comprehensive Genomic Report Dashboard</title>
        <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">
        <style>
            body {{ font-family: 'Arial', sans-serif; line-height: 1.6; padding: 20px; background-color: #f9f9f9; }}
            .container {{ background-color: white; padding: 30px; border-radius: 8px; box-shadow: 0 0 20px rgba(0,0,0,0.1); }}
            h1 {{ color: #0056b3; border-bottom: 3px solid #0056b3; padding-bottom: 10px; }}
            h2 {{ margin-top: 40px; color: #007bff; border-left: 5px solid #007bff; padding-left: 10px; }}
            h3 {{ color: #333; margin-top: 20px; }}
            .section-divider {{ margin: 50px 0; border-bottom: 2px dashed #ccc; }}
        </style>
    </head>
    <body>
        <div class="container">
            <header class="text-center mb-5">
                <h1>🧬 Comprehensive Bioinformatics Report Dashboard</h1>
                <p class="lead">Multi-modal analysis integrating CRISPR guides, sequence conservation, functional scoring, and antibiotic resistance profiling.</p>
            </header>

            {guides_html}

            <hr class="section-divider">
            {conservation_html}

            <hr class="section-divider">
            {evo2_html}

            <hr class="section-divider">
            {card_html}

            <footer class="text-center mt-5 py-3 border-top">
                <p>&copy; Bioinformatics Analysis Dashboard | Generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </footer>
        </div>
    </body>
    </html>
    """
    return html_dashboard

# --- 4. Execution and Output ---

if __name__ == "__main__":
    
    # Ensure the directories exist for clean execution environment checks (if running standalone)
    # In a real environment, these directories must exist.
    if not os.path.exists(GUIDES_DIR):
        os.makedirs(GUIDES_DIR)
    if not os.path.exists(REPORT_DIR):
        os.makedirs(REPORT_DIR)
    if not os.path.exists(EVO2_DIR):
        os.makedirs(EVO2_DIR)


    # Generate the final HTML report
    final_report_html = create_html_dashboard(
        df_guides, 
        df_conservation, 
        df_evo2, 
        df_card, 
        chart_conservation, 
        chart_evo2
    )

    # Save the report
    OUTPUT_FILENAME = "comprehensive_report.html"
    with open(OUTPUT_FILENAME, 'w') as f:
        f.write(final_report_html)

    print(f"\n==================================================")
    print(f"SUCCESS: The comprehensive report has been saved to {OUTPUT_FILENAME}")
    print(f"==================================================")
