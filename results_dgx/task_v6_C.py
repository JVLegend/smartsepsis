import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import base64
from io import BytesIO

# --- 1. Setup ---
# Define file paths
CONSERVATION_FILE = 'reports/conservation_analysis.tsv'
EVO2_FILE = 'reports/evo2_scoring/functional_scores.tsv'
OUTPUT_HTML = 'dashboard.html'

# --- 2. Data Loading ---
print("--- Loading Data ---")

# Load Conservation Analysis Data (TSV)
try:
    df_con = pd.read_csv(CONSERVATION_FILE, sep='\t')
    print(f'Successfully loaded {CONSERVATION_FILE}')
    print(f'Columns: {df_con.columns.tolist()}')
except FileNotFoundError:
    # Create dummy data if the file is missing for demonstration purposes
    print(f"Warning: {CONSERVATION_FILE} not found. Using dummy data.")
    data = {
        'variant': ['V1', 'V1', 'V2', 'V3', 'V3'],
        'family': ['A', 'A', 'B', 'B', 'C'],
        'accession': ['ACC1', 'ACC1', 'ACC2', 'ACC2', 'ACC3'],
        'match': ['exact', 'near', 'exact', 'none', 'exact'],
        'mismatches': [0, 1, 0, 2, 0],
        'position': [10, 12, 20, 21, 30],
        'strand': ['+', '+', '-', '-', '+']
    }
    df_con = pd.DataFrame(data)

# Load Evo2 Functional Scores Data (TSV)
try:
    df_evo2 = pd.read_csv(EVO2_FILE, sep='\t')
    print(f'Successfully loaded {EVO2_FILE}')
    print(f'Columns: {df_evo2.columns.tolist()}')
except FileNotFoundError:
    # Create dummy data if the file is missing for demonstration purposes
    print(f"Warning: {EVO2_FILE} not found. Using dummy data.")
    data = {
        'variant': ['V1', 'V2', 'V3', 'V4', 'V5'],
        'reference': ['R1', 'R2', 'R3', 'R4', 'R5'],
        'functional_distance': [0.1, 0.8, 0.3, 0.1, 0.9],
        'impact': [5, 10, 3, 1, 8],
        'confidence': [0.9, 0.8, 0.7, 0.9, 0.6],
        'resistance_maintained': [True, False, True, True, False],
        'gc_shift': [2, 1, 3, 2, 1],
        'kmer_distance': [1, 2, 1, 3, 2]
    }
    df_evo2 = pd.DataFrame(data)

# --- 3. Analysis and Visualization ---

# --- CHART 1: Conservation Rates per Gene/Family ---
def plot_conservation(df):
    # Since 'gene' column is not explicitly available in df_con, we group by 'family' 
    # as the best proxy for gene/grouping unit available.
    print("\nGenerating Conservation Rate Plot...")
    
    # Calculate the percentage of exact matches per family
    con_rates = df.groupby('family')['match'].value_counts(normalize=True).unstack(fill_value=0)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot the relative proportions of match types (exact, near, none)
    con_rates.plot(kind='bar', ax=ax)
    
    plt.title('Conservation Rate Distribution by Family (Match Type)')
    plt.xlabel('Family')
    plt.ylabel('Proportion (%)')
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Match Type', loc='best')
    plt.tight_layout()
    
    return fig

fig_con = plot_conservation(df_con)


# --- CHART 2: Evo2 Scores Distribution ---
def plot_evo2(df):
    print("Generating Evo2 Distribution Plot...")
    
    # We will visualize the distribution of the functional distance score
    plt.figure(figsize=(10, 6))
    sns.histplot(df['functional_distance'], kde=True, bins=15)
    plt.title('Distribution of Functional Distance Scores (Evo2)')
    plt.xlabel('Functional Distance')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.5)
    plt.tight_layout()
    
    return plt.gcf()

fig_evo2 = plot_evo2(df_evo2)


# --- CHART 3: Variant Coverage Proxy ---
def plot_coverage(df_con, df_evo2):
    print("Generating Variant Coverage Plot...")
    
    # We will use the count of unique records in each dataset as a proxy for coverage magnitude
    
    count_con = df_con[['variant']].drop_duplicates().shape[0]
    count_evo2 = df_evo2[['variant']].drop_duplicates().shape[0]
    
    data_coverage = pd.DataFrame({
        'Dataset': ['Conservation Analysis', 'Evo2 Scoring'],
        'Unique Variant Count': [count_con, count_evo2]
    })
    
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.barplot(x='Dataset', y='Unique Variant Count', data=data_coverage, palette='viridis')
    
    plt.title('Proxy Variant Coverage: Count of Unique Variants')
    plt.xlabel('')
    plt.ylabel('Number of Unique Variants')
    plt.ylim(0, max(count_con, count_evo2) * 1.2)
    plt.tight_layout()
    
    return fig

fig_coverage = plot_coverage(df_con, df_evo2)


# --- 4. Dashboard Generation ---
def create_dashboard(fig_con, fig_evo2, fig_coverage, output_path):
    print("\nCreating Dashboard...")
    
    # Use matplotlib's figures directly to capture them in base64 format
    # We must close existing figures before creating the main combined figure
    plt.close(fig_con)
    plt.close(fig_evo2)
    plt.close(fig_coverage)
    
    # Create the combined figure layout
    fig_dashboard = plt.figure(figsize=(16, 12))
    
    # Add subplots
    gs = fig_dashboard.add_gridspec(3, 1, height_ratios=[1, 1, 1])
    
    ax1 = fig_dashboard.add_subplot(gs[0])
    ax2 = fig_dashboard.add_subplot(gs[1])
    ax3 = fig_dashboard.add_subplot(gs[2])
    
    # Recreate/copy plotting elements onto the combined figure structure for display
    # Note: This is complex because the previous fig objects hold the state. 
    # For simplicity and robustness, we will save the figures temporarily and embed the images.
    
    # Alternative approach: Save each figure to a buffer and display the image in HTML.
    
    # Save figures to memory buffers
    buf_con = BytesIO()
    fig_con.savefig(buf_con, format='png')
    buf_con.seek(0)
    base64_con = base64.b64encode(buf_con.read()).decode('utf-8')

    buf_evo2 = BytesIO()
    fig_evo2.savefig(buf_evo2, format='png')
    buf_evo2.seek(0)
    base64_evo2 = base64.b64encode(buf_evo2.read()).decode('utf-8')
    
    buf_coverage = BytesIO()
    fig_coverage.savefig(buf_coverage, format='png')
    buf_coverage.seek(0)
    base64_coverage = base64.b64encode(buf_coverage.read()).decode('utf-8')

    # Build the HTML content
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Bioinformatics Functional Analysis Dashboard</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1 {{ color: #333; border-bottom: 2px solid #eee; padding-bottom: 10px; }}
            .chart-container {{ margin-bottom: 40px; padding: 20px; border: 1px solid #ddd; border-radius: 8px; }}
            h2 {{ color: #555; }}
            img {{ max-width: 100%; height: auto; display: block; }}
        </style>
    </head>
    <body>
        <h1>Comprehensive Functional Genomics Dashboard</h1>
        
        <div class="chart-container">
            <h2>1. Conservation Rates by Family (Match Analysis)</h2>
            <p>Shows the proportional distribution of sequence match types (exact, near, none) across different evolutionary families.</p>
            <img src="data:image/png;base64,{base64_con}" alt="Conservation Rates Chart">
        </div>
        
        <div class="chart-container">
            <h2>2. Evo2 Functional Distance Distribution</h2>
            <p>Displays the distribution of functional distance scores, indicating how far the variant functionally deviates from the reference.</p>
            <img src="data:image/png;base64,{base64_evo2}" alt="Evo2 Distribution Chart">
        </div>
        
        <div class="chart-container">
            <h2>3. Variant Coverage Proxy (Unique Variants)</h2>
            <p>Compares the number of unique variants processed in the two main datasets.</p>
            <img src="data:image/png;base64,{base64_coverage}" alt="Variant Coverage Chart">
        </div>

        <footer>
            <p style="text-align: center; margin-top: 30px; color: #777;">Dashboard generated by Bioinformatics Analysis Script.</p>
        </footer>
    </body>
    </html>
    """

    # Save the combined HTML content
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    print(f"\n✅ Successfully created dashboard at: {output_path}")

# --- 5. Execution ---
if __name__ == "__main__":
    # Generate figures (these figures hold the necessary plot data)
    fig_con = plot_conservation(df_con)
    fig_evo2 = plot_evo2(df_evo2)
    fig_coverage = plot_coverage(df_con, df_evo2)
    
    # Create and save the final HTML dashboard
    create_dashboard(fig_con, fig_evo2, fig_coverage, OUTPUT_HTML)
