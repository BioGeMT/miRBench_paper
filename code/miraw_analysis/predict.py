import pandas as pd
import numpy as np
from collections import Counter
from itertools import product
import pickle
import argparse
import logging
import os

def setup_logging(log_file):
    """
    Set up logging to write logs to a file and display them on the console.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),    # Write logs to the specified file
            logging.StreamHandler()          # Show logs in the console
        ]
    )

def get_all_possible_kmers(k):
    """Generate all possible k-mers for a given k."""
    return [''.join(kmer) for kmer in product("ACGT", repeat=k)]

def kmer_count_matrix(sequences, k):
    """Create a k-mer count matrix for a list of sequences."""
    all_kmers = get_all_possible_kmers(k)
    kmer_counts = []
    for sequence in sequences:
        kmer_count = Counter([sequence[i:i+k] for i in range(len(sequence) - k + 1)])
        row = [kmer_count.get(kmer, 0) for kmer in all_kmers]
        kmer_counts.append(row)
    return pd.DataFrame(kmer_counts, columns=all_kmers)

def generate_random_predictions(test_df):
    """
    Generates random predictions between 0 and 1 for each sample in X_test.
    
    Args:
    - X_test: The test set (features)
    
    Returns:
    - random_preds: Array of random predictions rounded to 4 decimal places
    """
    # Set the seed for reproducibility (fixed to 42)
    np.random.seed(42)
    
    # Generate random predictions between 0 and 1
    random_preds = np.random.rand(len(test_df))
    
    # Round the predictions to 4 decimal places
    random_preds = np.round(random_preds, 4)
    
    return random_preds

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Load the trained Decision Tree classifier and predict on test set.")
    parser.add_argument("--test_set", type=str, required=True, help="Path to the test dataset (.tsv).")
    parser.add_argument("--k", type=int, required=True, help="Length of k-mers to use for feature extraction.")
    parser.add_argument("--models", nargs='+', required=True, help="List of paths to the trained models (.pkl).")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the predictions (.tsv).")
    args = parser.parse_args()
    
    # Set up logging
    test_set_name = os.path.splitext(os.path.basename(args.test_set))[0]
    log_file = f"{args.output_dir}/{test_set_name}_predictions.log"
    setup_logging(log_file)

    # Load datasets
    logging.info(f"Loading test set from {args.test_set}...")
    test_df = pd.read_csv(args.test_set, sep="\t")

    # Check for miRAW random/regular Split column names and rename to standardise to miRBench datasets column names
    if "Mature_mirna_transcript" in test_df.columns and "Positive_Negative" in test_df.columns:
        test_df = test_df.rename(columns={"Mature_mirna_transcript": "noncodingRNA", "Positive_Negative": "label"})

    # Check for miRAW all Training Sites column names and rename to standardise to miRBench datasets column names
    if "mature_miRNA_Transcript" in test_df.columns and "validation" in test_df.columns:
        test_df = test_df.rename(columns={"mature_miRNA_Transcript": "noncodingRNA", "validation": "label"})

    # Generate k-mer count matrices
    logging.info(f"Generating k-mer features for k={args.k}...")
    X_test = kmer_count_matrix(test_df["noncodingRNA"], args.k)
    y_test = test_df["label"]

    for model_path in args.models:
        # Load the trained saved model
        logging.info(f"Loading the trained model from {model_path}...")
        with open(model_path, "rb") as model_file:
            model = pickle.load(model_file)

        # Predict on the test set
        logging.info("Predicting on the test set...")
        y_preds = model.predict_proba(X_test)[:, 1]

        # Append predictions to the test set
        model_name = os.path.splitext(os.path.basename(model_path))[0] + '_model'
        test_df[model_name] = y_preds

    # Generate random predictions
    logging.info("Generating predictions for a random model...")
    random_preds = generate_random_predictions(test_df)

    # Add the random predictions to the test set
    test_df['Random_model'] = random_preds
    
    # Save the predictions
    output_file_path = f"{args.output_dir}/{test_set_name}_predictions.tsv"
    logging.info(f"Saving all predictions to {output_file_path}...")
    test_df.to_csv(output_file_path, sep="\t", index=False)
    
    logging.info("Done!")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception("An error occurred during execution.")
        raise
