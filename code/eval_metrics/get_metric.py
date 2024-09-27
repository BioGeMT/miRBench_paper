import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
import argparse
import sys

def load_data(input_file):
    # load the data from the input file
    return pd.read_csv(input_file, sep='\t')

def get_metric(data, predictors, metric):
        
    metric_dict = {}

    for predictor in predictors:
        if predictor not in data.columns:
            raise KeyError(f"Predictor {predictor} not found in the data.")

        if metric == 'auc_pr':
            if predictor.startswith('Seed'):
                continue
            else:        
                precision, recall, _ = precision_recall_curve(data['label'], data[predictor])
                pr_auc = auc(recall, precision)
                metric_dict[predictor] = np.round(pr_auc, 2)

        elif metric == 'auc_roc':
            if predictor.startswith('Seed'):
                continue
            else:
                roc_auc = roc_auc_score(data['label'], data[predictor])
                metric_dict[predictor] = np.round(roc_auc, 2)

        elif metric == 'avg_p_score':
            if predictor.startswith('Seed'):
                continue
            else:
                precision, recall, _ = precision_recall_curve(data['label'], data[predictor])
                avg_p_score = average_precision_score(data['label'], data[predictor])
                metric_dict[predictor] = np.round(avg_p_score, 2)

        else:
            raise ValueError(f"Invalid metric: {metric}. Please choose one of 'auc-pr', 'auc-roc', or 'avg_p_score'.")

    return metric_dict

###############

def main():
    # argument parsing
    parser = argparse.ArgumentParser(description="Evaluate predictors.")
    parser.add_argument('--ifile', help="Input file containing the prediction scores in TSV format (default: STDIN)", default=None)
    parser.add_argument('--predictors', help="List of predictor names (default: all)", default=None)
    parser.add_argument('--metric', help="Evaluation metric to compute; auc_pr, auc_roc, or avg_p_score.", default=None)
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
        args.predictors = ['TargetScanCnn_McGeary2019',
            'CnnMirTarget_Zheng2020',
            'TargetNet_Min2021',
            'miRBind_Klimentova2022',
            'miRNA_CNN_Hejret2023',
            'InteractionAwareModel_Yang2024', 
            'RNACofold',
            'Seed8mer', 
            'Seed7mer', 
            'Seed6mer', 
            'Seed6merBulgeOrMismatch'] # In chronological order

    # if metric is none, raise an error
    if args.metric is None:
        raise ValueError(f"Missing metric. Please choose one of 'auc_pr', 'auc_roc', or 'avg_p_score'.")

    # get the metrics
    metric = get_metric(data, args.predictors, args.metric)

    # write the results to the output file
    with open(args.ofile, 'w') as ofile:
        ofile.write(f"Tool\t{args.metric}\n")
        for predictor in args.predictors:
            if predictor.startswith('Seed'):
                continue
            else:   
                ofile.write(f"{predictor}\t{metric[predictor]}\n")

if __name__ == "__main__":
    main()
