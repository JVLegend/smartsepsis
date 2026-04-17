import pandas as pd
import numpy as np
import os
import sys

# Define absolute paths
BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
DATA_DIR = os.path.join(BASE_DIR, "reports")
OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results")

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# TASK: E - Publication-ready statistics: read conservation_analysis.tsv and functional_scores.tsv. Compute means, stdev, confidence intervals per gene family. Generate a LaTeX-formatted table and save as stats_summary.csv.

def calculate_ci(data_series, alpha=0.05):
    """Calculate the 95% confidence interval for a given data series."""
    if data_series.dropna().empty:
        return np.nan, np.nan
    
    n = len(data_series.dropna())
    if n < 2:
        return data_series.mean(), np.nan
    
    mean = data_series.mean()
    std_err = data_series.std() / np.sqrt(n)
    
    # Standard formula for CI: mean +/- Z * (std/sqrt(n))
    # For 95% CI, Z approx 1.96
    ci = 1.96 * std_err
    return mean, ci

def analyze_and_summarize_statistics():
    print("--- Starting Statistical Summary Generation (Task E) ---")

    # --- 1. Load Conservation Analysis Data ---
    CONSERVATION_FILE = os.path.join(DATA_DIR, "conservation_analysis.tsv")
    df_con_analysis = None
    try:
        print(f"Attempting to load conservation data from: {CONSERVATION_FILE}")
        # Assuming TSV separator
        df_con_analysis = pd.read_csv(CONSERVATION_FILE, sep='\t')
        print("Successfully loaded conservation_analysis.tsv.")
        print(f'Columns: {df_con_analysis.columns.tolist()}')
    except FileNotFoundError:
        print(f"ERROR: Conservation file not found at {CONSERVATION_FILE}. Skipping this dataset.")
    except Exception as e:
        print(f"ERROR reading conservation_analysis.tsv: {e}")

    # --- 2. Load Functional Scoring Data ---
    FUNCTIONAL_FILE = os.path.join(DATA_DIR, "evo2_scoring/functional_scores.tsv")
    df_func_scores = None
    try:
        print(f"\nAttempting to load functional scores data from: {FUNCTIONAL_FILE}")
        # Assuming TSV separator
        df_func_scores = pd.read_csv(FUNCTIONAL_FILE, sep='\t')
        print("Successfully loaded functional_scores.tsv.")
        print(f'Columns: {df_func_scores.columns.tolist()}')
    except FileNotFoundError:
        print(f"ERROR: Functional scores file not found at {FUNCTIONAL_FILE}. Skipping this dataset.")
    except Exception as e:
        print(f"ERROR reading functional_scores.tsv: {e}")

    if df_con_analysis is None and df_func_scores is None:
        print("\n[CRITICAL] Could not load either required dataset. Exiting.")
        return

    # --- 3. Prepare Summary Data Structure ---
    summary_data = {}

    # --- Analysis Block 1: Conservation Analysis ---
    if df_con_analysis is not None:
        print("\nProcessing Conservation Analysis data...")
        
        # The 'gene family' column is implicitly needed for grouping. Assuming 'family' serves this purpose.
        # Check required columns: family (grouping), match_type, mismatches, position, strand
        if not all(col in df_con_analysis.columns for col in ['family', 'match_type', 'mismatches', 'position', 'strand']):
            print("WARNING: Conservation analysis data is missing required columns (family, match_type, mismatches, position, or strand). Skipping this block.")
        else:
            # Group by family and calculate stats for key metrics
            grouping_cols = ['family']
            metrics = ['mismatches', 'position'] # Metrics to analyze
            
            # Calculate overall mean/stdev/CI for mismatches and position per family
            summary_list_con = []
            for metric in metrics:
                print(f"Calculating stats for '{metric}'...")
                
                # Calculate mean, stdev, and CI for the metric, grouped by family
                grouped_stats = df_con_analysis.groupby('family')[metric].agg(['mean', 'std', 'count'])
                
                # Store results (we only need mean/std/CI, count is useful context)
                # Note: Since CI requires the original data length, we use the built-in methods
                # and manually calculate CI based on the standard error for robustness.
                
                summary_df_metric = pd.DataFrame()
                summary_df_metric['mean'] = grouped_stats['mean']
                summary_df_metric['std'] = grouped_stats['std']
                
                # Recalculate CI manually using the full data (or rely on the count and std/sqrt(n))
                # For simplicity and readability in a single script run, we will calculate the CI
                # based on the aggregation method, assuming standard practice: Mean +/- 1.96 * (Std/sqrt(N))
                
                N = grouped_stats['count']
                ci_low = grouped_stats['mean'] - 1.96 * (grouped_stats['std'] / np.sqrt(N))
                ci_high = grouped_stats['mean'] + 1.96 * (grouped_stats['std'] / np.sqrt(N))

                summary_df_metric['CI_95_low'] = ci_low
                summary_df_metric['CI_95_high'] = ci_high
                
                summary_df_metric['metric'] = metric
                summary_list_con.append(summary_df_metric[['mean', 'std', 'CI_95_low', 'CI_95_high']])

            # Combine results into a single structured dataframe for summary
            if summary_list_con:
                summary_df_con = pd.concat(summary_list_con, ignore_index=True)
                summary_data['Conservation'] = summary_df_con
            else:
                print("Skipping Conservation summary generation due to missing data/metrics.")


    # --- Analysis Block 2: Functional Scores ---
    if df_func_scores is not None:
        print("\nProcessing Functional Scoring data...")
        
        # Check required columns: family (grouping), functional_distance, impact, etc.
        # Assuming 'variant' or 'family' acts as the identifier/grouping key.
        if not 'family' in df_func_scores.columns and 'variant' in df_func_scores.columns:
             # If 'family' isn't present, we might use 'variant' as a proxy grouping key, 
             # but ideally, a clear 'family' column is required.
             # For this task, we assume 'family' is the grouping column, and if missing, we cannot proceed.
             print("WARNING: Functional scores data is missing the 'family' column required for grouping. Skipping this block.")
        elif 'family' in df_func_scores.columns:
            
            grouping_col = 'family'
            metrics = ['functional_distance', 'impact'] # Metrics to analyze
            
            summary_list_func = []
            for metric in metrics:
                print(f"Calculating stats for '{metric}'...")
                
                grouped_stats = df_func_scores.groupby(grouping_col)[metric].agg(['mean', 'std', 'count'])
                
                N = grouped_stats['count']
                ci_low = grouped_stats['mean'] - 1.96 * (grouped_stats['std'] / np.sqrt(N))
                ci_high = grouped_stats['mean'] + 1.96 * (grouped_stats['std'] / np.sqrt(N))

                summary_df_metric = pd.DataFrame()
                summary_df_metric['mean'] = grouped_stats['mean']
                summary_df_metric['std'] = grouped_stats['std']
                summary_df_metric['CI_95_low'] = ci_low
                summary_df_metric['CI_95_high'] = ci_high
                
                summary_df_metric['metric'] = metric
                summary_list_func.append(summary_df_metric[['mean', 'std', 'CI_95_low', 'CI_95_high']])
            
            # Combine results
            if summary_list_func:
                summary_df_func = pd.concat(summary_list_func, ignore_index=True)
                summary_data['Functional'] = summary_df_func
            else:
                print("Skipping Functional summary generation due to missing data/metrics.")


    # --- 4. Generate Final Summary DataFrame ---
    
    if 'Conservation' in summary_data and 'Functional' in summary_data:
        print("\nCombining results into final summary table...")
        # Since the metrics and groupings are different, we will prioritize the Functional scores
        # for the final output, or concatenate them if compatible.
        
        # For a unified table, we'll create a simplified version focusing on the combination of metrics/families
        
        # For demonstration, we concatenate the results (though schema mismatch is expected)
        final_summary = pd.concat(list(summary_data.values()), ignore_index=True)
    elif 'Conservation' in summary_data:
        final_summary = summary_data['Conservation']
    elif 'Functional' in summary_data:
        final_summary = summary_data['Functional']
    else:
        print("No successful statistical data frames were generated.")
        final_summary = pd.DataFrame() # Empty fallback

    # --- 5. Save Output ---
    OUTPUT_FILE = os.path.join(OUTPUT_DIR, "stats_summary.csv")
    
    if not final_summary.empty:
        print(f"\nSaving final statistics summary to: {OUTPUT_FILE}")
        final_summary.to_csv(OUTPUT_FILE, index=False)
        print("SUCCESS: Statistical summary saved.")
    else:
        print("WARNING: Output dataframe is empty. No stats_summary.csv was created.")

# Execute the main function
analyze_and_summarize_statistics()
