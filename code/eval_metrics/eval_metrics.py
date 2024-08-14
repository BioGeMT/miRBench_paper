import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import roc_auc_score
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
import argparse
import sys

def load_data(input_file):
    # load the data from the input file
    return pd.read_csv(input_file, sep='\t')

def get_metrics(data, predictors):

    f1score_dict = {}
    pr_auc_dict = {}
    roc_auc_dict = {}
    avg_p_score_dict = {}

    for predictor in predictors:
        if predictor not in data.columns:
            raise KeyError(f"Predictor {predictor} not found in the data.")
        if predictor.startswith('Seed'):
            _, _, f1score, _ = precision_recall_fscore_support(data['label'].values, data[predictor].values, average='binary', beta=1.0)
            f1score_dict[predictor] = f1score
        else:        
            precision, recall, _ = precision_recall_curve(data['label'], data[predictor])
            pr_auc = auc(recall, precision)
            pr_auc_dict[predictor] = pr_auc

            roc_auc = roc_auc_score(data['label'], data[predictor])
            roc_auc_dict[predictor] = roc_auc

            avg_p_score = average_precision_score(data['label'], data[predictor])
            avg_p_score_dict[predictor] = avg_p_score

    return pr_auc_dict, roc_auc_dict, avg_p_score_dict, f1score_dict

###############

def main():
    # argument parsing
    parser = argparse.ArgumentParser(description="Evaluate predictors using PR AUC, ROC AUC, average precision, F1 score")
    parser.add_argument('--ifile', help="Input file (default: STDIN)", default=None)
    parser.add_argument('--predictors', help="List of predictor names (default: all)", default=None)
    parser.add_argument('--ofile', help="Output file (default: STDOUT)", default=None)
    args = parser.parse_args()

    # if ifile is none, set it to sys.stdin
    if args.ifile is None:
        args.ifile = sys.stdin

    # load the data
    data = load_data(args.ifile)

    # if ofile is none, set it to sys.stdout
    if args.ofile is None:
        args.ofile = sys.stdout

    # if predictors is none, set it to all columns, except columns gene, noncodingRNA, noncodingRNA_fam, feature, and label
    if args.predictors is None:
        args.predictors = [col for col in data.columns if col not in ['gene', 'noncodingRNA', 'noncodingRNA_fam', 'feature', 'label']]

    # get the metrics
    pr_aucs, roc_aucs, avg_p_scores, fscore_dict = get_metrics(data, args.predictors)

    # write the results to the output file
    with open(args.ofile, 'w') as ofile:
        ofile.write("Predictor\tPR_AUC\tROC_AUC\tAverage_Precision\tF1_Score\n")
        for predictor in args.predictors:
            ofile.write(f"{predictor}\t{pr_aucs.get(predictor, 'NA')}\t{roc_aucs.get(predictor, 'NA')}\t{avg_p_scores.get(predictor, 'NA')}\t{fscore_dict.get(predictor, 'NA')}\n")

if __name__ == "__main__":
    main()
