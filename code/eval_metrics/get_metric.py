import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import confusion_matrix
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

        elif metric == 'cm':
            if predictor =='TargetScanCnn_McGeary2019' or predictor == 'RNACofold':
                # Normalise the predictions to [0, 1]
                scaler = MinMaxScaler()
                y_pred_reshaped = data[predictor].values.reshape(-1, 1)
                y_pred_normalised = scaler.fit_transform(y_pred_reshaped)
                y_pred = y_pred_normalised.flatten()
            else:
                y_pred = data[predictor].tolist()
            
            y_true = data['label'].tolist()

            if predictor.startswith('Seed'):
                y_pred_bin = y_pred
            else:
                # Compute the binary predictions
                precision, recall, thresholds = precision_recall_curve(y_true, y_pred) # some recall values are 0
                np.seterr(invalid='ignore') # ignore division by zero warning
                fscore = (2 * precision * recall) / (precision + recall)
                fscore_max_index = np.argmax(fscore) # locate the index of the largest f score
                threshold = thresholds[fscore_max_index]
                y_pred_bin = [1 if p >= threshold else 0 for p in y_pred]

            # Compute and extract TP, TN, FP, FN
            tn, fp, fn, tp = confusion_matrix(y_true, y_pred_bin).ravel()
            
            metric_dict[predictor] = [int(tn), int(fp), int(fn), int(tp)]

        else:
            raise ValueError(f"Invalid metric: {metric}. Please choose one of 'auc-pr', 'auc-roc', 'avg_p_score', or 'cm'.")

    return metric_dict

###############

def main():
    # argument parsing
    parser = argparse.ArgumentParser(description="Evaluate predictors.")
    parser.add_argument('--ifile', help="Input file containing the prediction scores in TSV format (default: STDIN)", default=None)
    parser.add_argument('--predictors', help="List of predictor names (default: all)", default=None)
    parser.add_argument('--metric', help="Evaluation metric to compute; auc_pr, auc_roc, avg_p_score, or cm.", default=None)
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
        raise ValueError(f"Missing metric. Please choose one of 'auc_pr', 'auc_roc', 'avg_p_score', or 'cm'.")

    # get the metrics
    metric = get_metric(data, args.predictors, args.metric)

    # write the results to the output file
    with open(args.ofile, 'w') as ofile:
        if args.metric == 'cm':
            ofile.write(f"Tool\tTN\tFP\tFN\tTP\n")
            for predictor in args.predictors:
                ofile.write(f"{predictor}\t{metric[predictor][0]}\t{metric[predictor][1]}\t{metric[predictor][2]}\t{metric[predictor][3]}\n")
        else:
            ofile.write(f"Tool\t{args.metric}\n")
            for predictor in args.predictors:
                if predictor.startswith('Seed'):
                    continue
                else:   
                    ofile.write(f"{predictor}\t{metric[predictor]}\n")

if __name__ == "__main__":
    main()
