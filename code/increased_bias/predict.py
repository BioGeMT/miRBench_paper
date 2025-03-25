import pandas as pd
import numpy as np
import argparse
import logging
import os

from miRBench.encoder import get_encoder
from miRBench.predictor import list_predictors, get_predictor

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

def generate_random_predictions(test_df):
    """
    Generates random predictions between 0 and 1 for each sample in X_test.
    
    Args:
    - X_test: The test set (features)
    
    Returns:
    - random_preds: Array of random predictions rounded to 4 decimal places
    """
    # Set the seed for reproducibility (fixed to 42)
    # np.random.seed(42)
    
    # Generate random predictions between 0 and 1
    random_preds = np.random.rand(len(test_df))
    
    # Round the predictions to 4 decimal places
    random_preds = np.round(random_preds, 4)
    
    return random_preds

def main():
    parser = argparse.ArgumentParser(description="Evaluate the miRBench models on the test set.")
    parser.add_argument("--test", required=True, help="Path to the test set.")
    parser.add_argument("--output", required=True, help="Path to the output file.")
    args = parser.parse_args()

    # Set up logging
    log_file = os.path.join(os.path.dirname(args.output), "predict.log")
    setup_logging(log_file)

    # Load the test set
    test = pd.read_csv(args.test, sep="\t")

    for tool in list_predictors():
        logging.info(f"Predicting labels using {tool}...")
        encoder = get_encoder(tool)
        predictor = get_predictor(tool)
        input = encoder(test)
        test[f"{tool}"] = predictor(input)

    test["Random"] = generate_random_predictions(test)

    test.to_csv(args.output, sep="\t", index=False)
    logging.info(f"Predictions saved to {args.output}")

if __name__ == "__main__":
    main()