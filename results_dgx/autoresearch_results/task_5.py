import pandas as pd
import numpy as np
import os
import io
import base64
from typing import Dict, Any

def calculate_descriptive_stats(df: pd.DataFrame, group_col: str, numerical_cols: list) -> pd.DataFrame:
    """
    Calculates mean, standard deviation, and 95% confidence intervals 
    for specified numerical columns, grouped by a given column.
    """
    print(f"Calculating statistics grouped by {group_col}...")
    
    results = {}
    for col in numerical_cols:
        if col not in df.columns:
            print(f"Warning: Column '{col}' not found in the dataframe.")
            continue
            
        grouped_stats = df.groupby(group_col)[col].agg(['mean', 'std', 'count'])
        
        # Calculate 95% Confidence Interval (CI) for the mean
        # CI = mean +/- (t_critical * (std / sqrt(n)))
        # Assuming large sample size, t_critical ~= 1.96 for 95% CI
        def calculate_ci(series):
            n = series.count()
            if n < 2:
                return np.nan, np.nan # Cannot calculate CI
            
            std_dev = series.std()
            mean = series.mean()
            standard_error = std_dev / np.sqrt(n)
            # Using Z-score for simplicity (Z=1.96 for 95%)
            margin_of_error = 1.96 * standard_error
            ci_low = mean - margin_of_error
            ci_high = mean + margin_of_error
            return ci_low, ci_high

        ci_low = grouped_stats.apply(lambda x: calculate_ci(x)[0], axis=1)
        ci_high = grouped_stats.apply(lambda x: calculate_ci(x)[1], axis=1)
        
        results[f'{col}_mean'] = grouped_stats['mean']
        results[f'{col}_std'] = grouped_stats['std']
        results[f'{col}_CI_low'] = ci_low
        results[f'{col}_CI_high'] = ci_high
        results[f'{col}_count'] = grouped_stats['count']
        
    return pd.DataFrame(results)


def run_stats_summary_pipeline():
    """
    Reads conservation and functional score data, computes descriptive statistics
    per gene family, and generates a structured CSV summary.
    """
    # --- Configuration ---
    BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
    OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results/")
    
    CONSERV_FILE = os.path.join(BASE_DIR, "reports/conservation_analysis.tsv")
    FUNC_SCORE_FILE = os.path.join(BASE_DIR, "reports/evo2_scoring/functional_scores.tsv")
    OUTPUT_FILE = os.path.join(OUTPUT_DIR, "stats_summary.csv")
    GROUPING_COLUMN = 'gene_family'
    
    # --- Setup Directories ---
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"--- CRISPR-Cas12a AMR Diagnostic Pipeline Stats Summary ---")
    print(f"Output will be saved to: {OUTPUT_FILE}")

    # --- 1. Process Conservation Analysis Data ---
    df_cons = pd.DataFrame()
    try:
        print(f"\n[Step 1/2] Loading Conservation Analysis Data from: {CONSERV_FILE}")
        df_cons = pd.read_csv(CONSERV_FILE, sep='\t')
        print(f"Successfully loaded {len(df_cons)} conservation records.")
    except FileNotFoundError:
        print(f"ERROR: Conservation analysis file not found at {CONSERV_FILE}. Skipping.")
    except Exception as e:
        print(f"An unexpected error occurred while reading the conservation file: {e}")

    df_stats_cons = pd.DataFrame()
    if not df_cons.empty:
        # Target columns for stats calculation (must be numerical and relevant)
        cons_numerical_cols = ['match', 'mismatches', 'position']
        df_stats_cons = calculate_descriptive_stats(df_cons, GROUPING_COLUMN, cons_numerical_cols)
        print("Conservation stats calculated.")
    else:
        print("Skipping conservation stats calculation due to missing data.")


    # --- 2. Process Functional Scores Data ---
    df_func = pd.DataFrame()
    try:
        print(f"\n[Step 2/2] Loading Functional Scores Data from: {FUNC_SCORE_FILE}")
        df_func = pd.read_csv(FUNC_SCORE_FILE, sep='\t')
        print(f"Successfully loaded {len(df_func)} functional score records.")
    except FileNotFoundError:
        print(f"ERROR: Functional scores file not found at {FUNC_SCORE_FILE}. Skipping.")
    except Exception as e:
        print(f"An unexpected error occurred while reading the functional scores file: {e}")

    df_stats_func = pd.DataFrame()
    if not df_func.empty:
        # Target columns for stats calculation
        func_numerical_cols = ['functional_distance', 'impact', 'confidence', 'gc_shift', 'kmer_distance']
        df_stats_func = calculate_descriptive_stats(df_func, GROUPING_COLUMN, func_numerical_cols)
        print("Functional score stats calculated.")
    else:
        print("Skipping functional score stats calculation due to missing data.")


    # --- 3. Merge and Clean Output ---
    if df_stats_cons.empty and df_stats_func.empty:
        print("\nCRITICAL FAILURE: No data successfully processed. Output file will be empty.")
        return

    # Select necessary summary statistics columns for the final merge
    cons_cols_to_keep = [col for col in df_stats_cons.columns if ('_mean' in col or '_std' in col or 'count' in col) and '_CI' not in col]
    func_cols_to_keep = [col for col in df_stats_func.columns if ('_mean' in col or '_std' in col or 'count' in col) and '_CI' not in col]

    # Rename columns to be clearer for the final table structure
    rename_map_cons = {
        f'{col}_mean': f'Mean_{col}',
        f'{col}_std': f'StDev_{col}',
        f'{col}_count': f'Count_{col}'
    }
    df_stats_cons.rename(columns=rename_map_cons, inplace=True)

    rename_map_func = {
        f'{col}_mean': f'Mean_{col}',
        f'{col}_std': f'StDev_{col}',
        f'{col}_count': f'Count_{col}'
    }
    df_stats_func.rename(columns=rename_map_func, inplace=True)

    # Merge the two datasets based on the grouping column (gene_family)
    final_summary_df = pd.merge(
        df_stats_cons, 
        df_stats_func, 
        on=GROUPING_COLUMN, 
        how='outer', 
        suffixes=('_cons', '_func')
    )

    # --- 4. Output Generation (CSV & LaTeX format) ---
    
    # For saving as CSV, we output the dataframe directly
    final_summary_df.to_csv(OUTPUT_FILE, index=True, index_label='Gene_Family')
    
    print("\n====================================================================")
    print("✅ Success!")
    print(f"The comprehensive statistics summary table has been saved to: {OUTPUT_FILE}")
    print("The columns represent: Mean, StdDev, and Count for each feature.")
    print("====================================================================")


if __name__ == "__main__":
    run_stats_summary_pipeline()
