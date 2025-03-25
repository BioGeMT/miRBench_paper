import pandas as pd
import numpy as np
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

def make_biased_dataset(train, test):
    """
    Create a biased test dataset.
    Keep only those positive samples from the test set that miRNA from these samples are more frequent in the positive class in the training set.
    Keep only those negative samples from the test set that miRNA from these samples are more frequent in the negative class in the training set.
    """

    # Calculate miRNA counts for each class
    miRNA_counts_class = train.groupby(['noncodingRNA', 'label']).size().unstack()

    # go over samples in test set and assign label based on majority class in train set
    def predict_label(row):
        miRNA = row['noncodingRNA']
        if miRNA in miRNA_counts_class.index:
            if miRNA_counts_class.loc[miRNA, 0] > miRNA_counts_class.loc[miRNA, 1]:
                return 0
            else:
                return 1
        else:
            return 0
        
    test['predicted_label'] = test.apply(predict_label, axis=1)

    # filter test to keep only samples where label equals predicted_label and 
    test_filtered = test[test['label'] == test['predicted_label']]
    test_filtered = test_filtered.drop(columns=['predicted_label'])
    test_filtered = test_filtered.reset_index(drop=True)

    # balance positive and negative samples
    test_filtered_pos = test_filtered[test_filtered['label'] == 1]
    test_filtered_neg = test_filtered[test_filtered['label'] == 0]
    max_samples = min(len(test_filtered_pos), len(test_filtered_neg))
    test_filtered = pd.concat([test_filtered_pos.head(max_samples), test_filtered_neg.head(max_samples)])
    test_filtered = test_filtered.reset_index(drop=True)

    return test_filtered

def main():
    parser = argparse.ArgumentParser(description="Create a biased dataset by increasing the number of samples of a specific class.")
    parser.add_argument("--input_train_file", type=str, help="Path to the input train TSV file.")
    parser.add_argument("--input_test_file", type=str, help="Path to the input test TSV file.")
    parser.add_argument("--output_file", type=str, help="Path to the output TSV file.")
    args = parser.parse_args()

    # Set up logging
    log_file = os.path.join(os.path.dirname(args.output_file), "make_biased_dataset.log")
    setup_logging(log_file)

    # Load the input CSV file
    logging.info(f"Loading the input train CSV file from {args.input_train_file}")
    train = pd.read_csv(args.input_train_file, sep="\t")
    logging.info(f"Loading the input test CSV file from {args.input_test_file}")
    test = pd.read_csv(args.input_test_file, sep="\t")

    biased_test = make_biased_dataset(train, test)
    logging.info(f"Biased dataset shape: {biased_test.shape}")

    # Save the biased dataset to a CSV file
    logging.info(f"Saving the biased dataset to {args.output_file}")
    biased_test.to_csv(args.output_file, index=False, sep="\t")

if __name__ == "__main__":
    main()