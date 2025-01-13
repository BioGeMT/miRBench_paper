import pandas as pd
from sklearn.metrics import precision_recall_curve, auc, average_precision_score
import matplotlib.pyplot as plt
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

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Load the predictions and calculate evaluation metrics.")
    parser.add_argument("--predictions", type=str, required=True, help="Path to the predictions file (.tsv).")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to output PR curves and evaluation metrics.")
    args = parser.parse_args()
    
    # Set up logging
    predictions_file_name = os.path.splitext(os.path.basename(args.predictions))[0]
    log_file = f"{args.output_dir}/{predictions_file_name}_evaluation.log"
    setup_logging(log_file)

    logging.info(f"Loading predictions from {args.predictions}...")
    predictions = pd.read_csv(args.predictions, sep="\t")

    # Check for miRAW random/regular Split column names and rename to standardise to miRBench datasets column names
    if "Mature_mirna_transcript" in predictions.columns and "Positive_Negative" in predictions.columns:
        predictions = predictions.rename(columns={"Mature_mirna_transcript": "noncodingRNA", "Positive_Negative": "label"})

    # Check for miRAW all Training Sites column names and rename to standardise to miRBench datasets column names
    if "mature_miRNA_Transcript" in predictions.columns and "validation" in predictions.columns:
        predictions = predictions.rename(columns={"mature_miRNA_Transcript": "noncodingRNA", "validation": "label"})

    y_test = predictions["label"]
    models = [col for col in predictions.columns if '_model' in col]
    av_prec_score_dict = {}
    pr_auc_dict = {}
    for model in models:
        logging.info(f"Evaluating predictions for {model}...")
        y_preds = predictions[model]
        av_prec_score = average_precision_score(y_test, y_preds)
        av_prec_score_dict[model] = av_prec_score
        precision, recall, _ = precision_recall_curve(y_test, y_preds)
        pr_auc = auc(recall, precision)
        pr_auc_dict[model] = pr_auc
    
        # Save the PR curve plot
        pr_curve_output_path = f"{args.output_dir}/{predictions_file_name}_{model}_pr_curve.png"
        logging.info(f"Saving the Precision-Recall Curve plot to {pr_curve_output_path}...")
        plt.figure(figsize=(8, 6))
        plt.plot(recall, precision, label=f"PR Curve (AUC = {pr_auc:.4f})")
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title(f"PRC for {model} on {predictions_file_name}")
        plt.legend(loc="lower left")
        plt.savefig(pr_curve_output_path)
    
    test_set_filename = predictions_file_name.removesuffix("_predictions")
    metrics = [av_prec_score_dict, pr_auc_dict]
    metric_names = ["avg_prec_score", "auc_pr"]
    for metric, metric_name in zip(metrics, metric_names):
        metric_df = pd.DataFrame({"Test_set": [test_set_filename]})
        for model, value in metric.items():
            metric_df[model] = [round(value, 4)]
        metric_output_path = f"{args.output_dir}/{metric_name}.tsv"
        write_header = not os.path.exists(metric_output_path)
        logging.info(f"Saving {metric_name} to {metric_output_path} (append mode: {not write_header})...")
        metric_df.to_csv(metric_output_path, sep="\t", index=False, mode='a', header=write_header)

    logging.info("Done!")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception("An error occurred during execution.")
        raise
    