import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import base64
from io import BytesIO

# --- Configuration and Paths ---
BASE_DIR = "/home/oftalmousp/jv-teste/harvard_bacend/"
SEQ_DIR = os.path.join(BASE_DIR, "sequences/")
OUTPUT_DIR = os.path.join(BASE_DIR, "autoresearch_results/")

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)


def read_fasta_sequences(directory):
    """
    Reads all FASTA files from the specified directory and returns a dictionary
    mapping sequence names to sequences.
    """
    print(f"Starting FASTA reading from: {directory}")
    sequences = {}
    file_paths = []
    
    # List all FASTA files
    try:
        file_paths = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".fasta")]
    except FileNotFoundError:
        print(f"Error: Directory {directory} not found.")
        return sequences
    
    if not file_paths:
        print("No FASTA files found in the directory.")
        return sequences

    for i, file_path in enumerate(file_paths):
        try:
            print(f"Processing file {i+1}/{len(file_paths)}: {os.path.basename(file_path)}")
            sequence_data = ""
            header = ""
            
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Start of a new sequence
                        if header:
                            sequences[header] = sequence_data
                        header = line[1:]  # Remove '>'
                            
                        # Initialize sequence data for the new header
                        sequence_data = ""
                    elif header:
                        # Append sequence line (assuming DNA bases)
                        sequence_data += line.upper()

            # Store the last sequence read
            if header:
                sequences[header] = sequence_data

        except Exception as e:
            print(f"Warning: Could not process file {file_path}. Error: {e}")
            continue
    
    print(f"Finished reading sequences. Total sequences loaded: {len(sequences)}")
    return sequences

def hamming_distance(seq1, seq2):
    """Calculates the Hamming distance between two sequences."""
    if len(seq1) != len(seq2):
        # In a real-world scenario, one might pad or truncate.
        # Here, we assume the sequences should be of equal length for accurate Hamming calculation.
        # If lengths differ, the distance is undefined/invalid for standard Hamming.
        return np.nan
    
    return sum(1 for i in range(len(seq1)) if seq1[i] != seq2[i])

def compute_distance_matrix(sequences):
    """
    Computes the pairwise Hamming distance for all sequences.
    Returns a DataFrame representing the distance matrix.
    """
    seq_names = list(sequences.keys())
    N = len(seq_names)
    print(f"\n--- Computing {N}x{N} Hamming Distance Matrix ---")
    
    distance_matrix = []
    for i in range(N):
        row = []
        for j in range(N):
            name1 = seq_names[i]
            name2 = seq_names[j]
            seq1 = sequences[name1]
            seq2 = sequences[name2]
            
            # Only calculate distance if lengths match
            if len(seq1) != len(seq2):
                 distance = np.nan
            else:
                distance = hamming_distance(seq1, seq2)
            
            row.append(distance)
        distance_matrix.append(row)
        
    # Create a DataFrame with proper column and index labels
    distance_df = pd.DataFrame(distance_matrix, index=seq_names, columns=seq_names)
    return distance_df

def save_distance_matrix(df, output_path):
    """Saves the DataFrame as a CSV file."""
    try:
        df.to_csv(output_path, index=True, header=True)
        print(f"\n[SUCCESS] Distance matrix saved successfully to: {output_path}")
    except Exception as e:
        print(f"[ERROR] Failed to save distance matrix: {e}")

def create_and_save_dendrogram(distance_df, output_html_path):
    """
    Performs hierarchical clustering and saves the result as an interactive HTML/Base64 image.
    """
    print("\n--- Generating Phylogenetic Dendrogram ---")
    
    # 1. Convert the square distance matrix into condensed form required by pdist
    # We only need the upper triangle (excluding the diagonal)
    distance_array = np.array(distance_df.fillna(np.nan).values)
    
    # Since pdist expects a condensed matrix (only unique pairs, row-major order)
    # We use the actual distance values for linkage.
    # Ensure we exclude NaN values if the matrix contained unequal length sequences.
    
    # Fill NaNs with a large arbitrary number for distance calculation 
    # (or better, only process the non-NaN subset, but linking requires a consistent matrix)
    # For robust execution, we will treat the distances based on the actual non-NaN pairs.
    
    # Filter out NaN columns/rows first if they exist
    valid_indices = np.where(~np.isnan(distance_array) * 1)
    valid_distance_array = distance_array[valid_indices[0], :][valid_indices[1]]
    
    if len(valid_distance_array) < 2:
         print("[WARNING] Not enough valid pairs of sequences to generate a dendrogram.")
         return
         
    # 2. Calculate linkage (clustering)
    # 'ward' is common, but 'average' or 'complete' might be better depending on biology. We use 'average'.
    Z = linkage(valid_distance_array, method='average')
    
    # 3. Generate the dendrogram structure
    plt.figure(figsize=(15, 10))
    dendrogram(Z, labels=distance_df.index.tolist(), leaf_rotation=90, leaf_points=False)
    plt.title('AMR Sequence Phylogenetic Dendrogram (Hamming Distance)')
    plt.xlabel('Sequences')
    plt.ylabel('Distance')
    plt.tight_layout()

    # 4. Save the plot to an in-memory buffer (BytesIO)
    buf = BytesIO()
    plt.savefig(buf, format='png')
    plt.close() # Close the plot to free memory

    # 5. Encode the image to Base64
    image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')

    # 6. Save the HTML file containing the base64 image
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Phylogenetic Dendrogram</title>
        <style>
            body {{ font-family: sans-serif; margin: 20px; }}
            img {{ max-width: 100%; height: auto; border: 1px solid #ccc; padding: 10px; }}
            h1 {{ color: #333; }}
        </style>
    </head>
    <body>
        <h1>CRISPR-Cas12a AMR Sequence Dendrogram</h1>
        <p>This dendrogram visualizes the genetic relatedness based on pairwise Hamming distances between the input sequences.</p>
        <img src="data:image/png;base64,{image_base64}" alt="Dendrogram Plot">
    </body>
    </html>
    """

    output_path = os.path.join(OUTPUT_DIR, "dendrogram_heatmap.html")
    try:
        with open(output_path, "w") as f:
            f.write(html_content)
        print(f"\n[SUCCESS] Dendrogram HTML saved successfully to: {output_path}")
    except Exception as e:
        print(f"[ERROR] Failed to write HTML file: {e}")


# =================================================================
# MAIN EXECUTION START
# =================================================================

# TASK: G - Phylogenetic distances: read all FASTA from sequences/, compute pairwise Hamming distance between sequences. Create a distance matrix CSV and an HTML dendrogram using scipy.cluster.hierarchy with matplotlib base64.

def main():
    print("=================================================================")
    print("TASK: G - CRISPR-Cas12a AMR Diagnostic Pipeline - Phylogenetic Distance Calculation")
    print("=================================================================")

    # 1. Read all sequences
    sequences_dict = read_fasta_sequences(SEQ_DIR)

    if not sequences_dict:
        print("Cannot proceed with distance calculation due to no valid sequences.")
        return

    # 2. Compute distance matrix
    distance_df = compute_distance_matrix(sequences_dict)
    
    # 3. Save the distance matrix
    distance_output_path = os.path.join(OUTPUT_DIR, "hamming_distance_matrix.csv")
    save_distance_matrix(distance_df, distance_output_path)
    
    # 4. Generate and save the dendrogram
    dendrogram_html_path = os.path.join(OUTPUT_DIR, "dendrogram_heatmap.html")
    create_and_save_dendrogram(distance_df, dendrogram_html_path)
    
    print("\n=================================================================")
    print("Pipeline execution complete. Results saved in:")
    print(OUTPUT_DIR)
    print("=================================================================")

if __name__ == "__main__":
    # Check for necessary libraries installation (Optional, but good practice)
    try:
        import scipy
        import matplotlib
        import pandas
        print("Required libraries found: numpy, pandas, matplotlib, scipy.")
    except ImportError:
        print("="*80)
        print("CRITICAL ERROR: Please ensure all necessary libraries are installed.")
        print("Run: pip install numpy pandas scipy matplotlib")
        print("="*80)
    else:
        main()
