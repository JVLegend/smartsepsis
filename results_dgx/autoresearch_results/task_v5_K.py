import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import base64
from io import BytesIO
from datetime import datetime

# --- Constants and Paths ---
WORKING_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
REPORT_DIR = os.path.join(WORKING_DIR, "reports")
COVARIANCE_DIR = os.path.join(REPORT_DIR, "covariance_probes")
OUTPUT_DIR = os.path.join(WORKING_DIR, "autoresearch_results")
OUTPUT_REPORT_PATH = os.path.join(OUTPUT_DIR, "comprehensive_report.html")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Helper Functions ---

def load_file_safely(file_path, sep=','):
    """Loads a file using try/except for robust handling."""
    print(f"Attempting to load file: {file_path}")
    try:
        if sep == 'tab':
            df = pd.read_csv(file_path, sep='\t')
        else: # Default to comma separator for CSV
            df = pd.read_csv(file_path, sep=sep)
        print(f"Successfully loaded data from {os.path.basename(file_path)}. Columns: {df.columns.tolist()}")
        return df
    except FileNotFoundError:
        print(f"WARNING: File not found at {file_path}. Returning None DataFrame.")
        return None
    except pd.errors.EmptyDataError:
        print(f"WARNING: File {os.path.basename(file_path)} is empty. Returning None DataFrame.")
        return None
    except Exception as e:
        print(f"ERROR loading {os.path.basename(file_path)}: {e}. Returning None DataFrame.")
        return None

def create_and_encode_plot(data, plot_type="hist", title="Distribution Plot"):
    """Generates a matplotlib plot, saves it to a BytesIO buffer, and returns a base64 string."""
    if data.empty:
        print("Skipping plot generation: Input data is empty.")
        return None

    # Set consistent style
    plt.style.use('ggplot')

    fig, ax = plt.subplots(figsize=(10, 6))
    
    try:
        if plot_type == "hist":
            # Assuming 'match_type' or a similar categorical column for histogram/count
            column = data.columns[0] if len(data.columns) > 0 else None
            if column:
                sns.countplot(y=data[column], order=data[column].value_counts().index, ax=ax)
                ax.set_title(f"{title} (Counts by {column})")
                ax.set_xlabel("Count")
                ax.set_ylabel(column)
            else:
                ax.text(0.5, 0.5, "No suitable column found for plotting.", ha='center')
        
        elif plot_type == "bar_score":
            # Assuming 'functional_distance' or 'impact' for bar chart
            column = data.columns[np.where(data.iloc[:, 0] == 'functional_distance')[0] if 'functional_distance' in data.columns else 0][0]
            sns.barplot(x=data[data.columns[0]], y=data.columns[1], ax=ax)
            ax.set_title(f"{title} (Functional Distance vs {data.columns[1]})")
            ax.set_xlabel(data.columns[0])
            ax.set_ylabel(data.columns[1])
        else:
            ax.text(0.5, 0.5, "Unsupported plot type.", ha='center')

        plt.tight_layout()
        
        # Save plot to buffer
        img_buffer = BytesIO()
        plt.savefig(img_buffer, format='png')
        plt.close(fig)
        
        base64_encoded_plot = base64.b64encode(img_buffer.getvalue()).decode('utf-8')
        return base64_encoded_plot

    except Exception as e:
        print(f"Error generating plot: {e}")
        return None

# --- Data Loading and Processing ---

print("\n--- Starting Data Loading Process ---")

# 1. Conservation Analysis
conservation_file = os.path.join(REPORT_DIR, "conservation_analysis.tsv")
df_conservation = load_file_safely(conservation_file, sep='tab')

# 2. Evo2 Scoring
evo2_file = os.path.join(REPORT_DIR, "evo2_scoring/functional_scores.tsv")
df_evo2 = load_file_safely(evo2_file, sep='\t')

# 3. CARD Analysis
card_file = os.path.join(REPORT_DIR, "reports/card_new_variants.csv")
df_card = load_file_safely(card_file, sep=',')

# 4. Covariance Probes
df_covariance_list = []
print("\n--- Processing Covariance Probes Directory ---")
for filename in os.listdir(COVARIANCE_DIR):
    if filename.endswith(".tsv"):
        filepath = os.path.join(COVARIANCE_DIR, filename)
        print(f"Loading covariance file: {filename}")
        df_cov = load_file_safely(filepath, sep='\t')
        if df_cov is not None:
            df_covariance_list.append(df_cov)

df_covariance = None
if df_covariance_list:
    df_covariance = pd.concat(df_covariance_list, ignore_index=True)
    print(f"Successfully combined {len(df_covariance_list)} covariance files. Total records: {len(df_covariance)}")


# --- Visualization Generation ---

print("\n--- Generating Visualizations ---")

# 1. Conservation Plot (Match Type Distribution)
conservation_plot_data = df_conservation.dropna(subset=['match_type'])
conservation_plot_base64 = create_and_encode_plot(
    conservation_plot_data, 
    plot_type="hist", 
    title="Conservation Analysis: Distribution of Match Types"
)

# 2. Evo2 Plot (Functional Distance Distribution)
# We select functional_distance and impact for visualization
evo2_plot_data = df_evo2[['functional_distance', 'impact']].dropna(subset=['functional_distance'])
evo2_plot_base64 = create_and_encode_plot(
    evo2_plot_data, 
    plot_type="bar_score", 
    title="Evo2 Scoring: Functional Distance vs Impact"
)

# 3. CARD Plot (Simple check: unique drug classes)
# We don't plot CARD, but we summarize unique classes for the report text.
unique_card_classes = df_card['drug_classes'].dropna().unique()

# 4. Covariance Plot (Example: Correlation check or simple count)
covariance_plot_base64 = None
if df_covariance is not None and len(df_covariance) > 0:
    # Just plotting the head data or a dummy plot since specific correlation columns are unknown
    cov_plot_data = df_covariance.head(10) 
    covariance_plot_base64 = create_and_encode_plot(
        cov_plot_data, 
        plot_type="hist", 
        title="Covariance Probes: Initial Feature Distribution Check"
    )


# --- HTML Report Generation ---

print("\n--- Generating HTML Report ---")

def create_html_section(title, data_df, plot_base64, description=None):
    """Generates the HTML content block for a specific report section."""
    html = f'<section id="{title.lower().replace(" ", "-")}">'
    html += f'<h2>{title}</h2>'
    
    if description:
        html += f'<p class="summary-text">{description}</p>'

    # Display plot
    if plot_base64:
        html += f'<h3>Visual Summary</h3>'
        html += f'<img src="data:image/png;base64,{plot_base64}" alt="{title} Plot" style="max-width: 100%; height: auto; border: 1px solid #ddd; padding: 10px;">'
    else:
        html += '<p>No visualization data available or required.</p>'
        
    # Display data preview
    if data_df is not None and not data_df.empty:
        html += '<h3>Data Snapshot (First 5 Rows)</h3>'
        # Converting DF to HTML table, replacing potential pandas/IPython HTML rendering
        html += data_df.head(5).to_html(classes='table table-striped', index=False)
    elif data_df is not None and data_df.empty:
        html += '<p>The dataset loaded, but it is empty.</p>'
    else:
        html += '<p>Data loading failed or no relevant data available.</p>'
        
    html += '</section><hr>'
    return html

# --- Report Compilation ---

html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CRISPR-Cas12a AMR Diagnostic Pipeline Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f4f7f9; color: #333; }}
        header {{ background-color: #007bff; color: white; padding: 20px; text-align: center; border-radius: 8px; margin-bottom: 30px; }}
        h1, h2, h3 {{ color: #0056b3; }}
        section {{ background-color: white; padding: 25px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .summary-text {{ background-color: #e9f5ff; padding: 15px; border-left: 5px solid #007bff; margin-top: 15px; font-style: italic; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 20px; }}
        th, td {{ padding: 10px 15px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #f2f2f2; color: #333; }}
        img {{ max-width: 100%; height: auto; display: block; margin: 20px 0; }}
    </style>
</head>
<body>

<header>
    <h1>🧬 CRISPR-Cas12a AMR Diagnostic Pipeline Report</h1>
    <p>Comprehensive Analysis of Guide Design, Resistance, and Evolutionary Markers</p>
    <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
</header>

"""

# 1. Conservation Section
description_conservation = "The conservation analysis tracks the relationship between the tested gene variant and its known counterpart. Results are categorized by match type (exact, near, none) to assess the conservation level of the guide region."
html_content += create_html_section(
    "I. Conservation Analysis (Variant Alignment)",
    df_conservation,
    conservation_plot_base64,
    description_conservation
)

# 2. Evo2 Scoring Section
description_evo2 = "Evo2 scoring provides quantitative metrics on the evolutionary stability and functional distance of the guide. High functional distance or low impact may indicate potential weakness or non-conserved regions."
html_content += create_html_section(
    "II. Evolutionary Scoring (Evo2)",
    df_evo2,
    evo2_plot_base64,
    description_evo2
)

# 3. CARD Section
card_description = (
    f"The CARD analysis identifies clinically relevant antibiotic resistance mechanisms. Based on the processed data, the variants are associated with multiple drug classes, including {', '.join(unique_card_classes[:5])}{'...' if len(unique_card_classes) > 5 else ''}."
)
html_content += f"""
<section id="card-analysis">
    <h2>III. Antibiotic Resistance Profile (CARD)</h2>
    <p class="summary-text">{card_description}</p>
    <h3>Data Snapshot (First 5 Rows)</h3>
    {df_card.head(5).to_html(classes='table table-striped', index=False)}
</section><hr>
"""

# 4. Covariance Probes Section
covariance_description = "Covariance probes provide biophysical features related to DNA structure (e.g., kmer composition, local melting temperature estimates). The aggregated data allows assessment of physical constraints on guide stability."
html_content += create_html_section(
    "IV. Biophysical Features (Covariance Probes)",
    df_covariance,
    covariance_plot_base64,
    covariance_description
)

# Final Save
with open(OUTPUT_REPORT_PATH, 'w') as f:
    f.write(html_content)

print("\n✅ Task K Complete. Comprehensive report saved successfully to:")
print(f"{OUTPUT_REPORT_PATH}")
