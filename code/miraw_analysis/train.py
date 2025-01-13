import pandas as pd
from collections import Counter
from itertools import product
from sklearn.tree import DecisionTreeClassifier
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

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Train a Decision Tree classifier using k-mer features and save the model.")
    parser.add_argument("--train_set", type=str, required=True, help="Path to the training dataset (.tsv).")
    parser.add_argument("--k", type=int, required=True, help="Length of k-mers.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the trained model (.pkl).")
    args = parser.parse_args()
    
    # Set up logging
    train_set_name = os.path.splitext(os.path.basename(args.train_set))[0]
    log_file = f"{args.output_dir}/DecisionTree_{args.k}mer_{train_set_name}.log"
    setup_logging(log_file)

    # Load datasets
    logging.info(f"Loading dataset from {args.train_set}...")
    train_df = pd.read_csv(args.train_set, sep="\t")

    # Check for miRAW random/regular Split column names and rename to standardise to miRBench datasets column names
    if "Mature_mirna_transcript" in train_df.columns and "Positive_Negative" in train_df.columns:
        train_df = train_df.rename(columns={"Mature_mirna_transcript": "noncodingRNA", "Positive_Negative": "label"})

    # Check for miRAW all Training Sites column names and rename to standardise to miRBench datasets column names
    if "mature_miRNA_Transcript" in train_df.columns and "validation" in train_df.columns:
        train_df = train_df.rename(columns={"mature_miRNA_Transcript": "noncodingRNA", "validation": "label"})

    # Generate k-mer count matrices
    logging.info(f"Generating k-mer features for k={args.k}...")
    X_train = kmer_count_matrix(train_df["noncodingRNA"], args.k)
    y_train = train_df["label"]
    
    # Train the Decision Tree model
    logging.info("Training the Decision Tree classifier...")
    model = DecisionTreeClassifier(random_state=42)
    model.fit(X_train, y_train)
    
    # Save the trained model
    output_file_path = f"{args.output_dir}/DecisionTree_{args.k}mer_{train_set_name}.pkl"
    logging.info(f"Saving the trained model to {output_file_path}...")
    with open(output_file_path, "wb") as f:
        pickle.dump(model, f)

    logging.info("Done!")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception("An error occurred during execution.")
        raise
