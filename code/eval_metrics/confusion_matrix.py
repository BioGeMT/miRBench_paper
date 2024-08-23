import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import confusion_matrix
import argparse
import sys

def load_data(input_file):
    # load the data from the input file
    return pd.read_csv(input_file, sep='\t')

def get_confusion_matrix(data, predictors):

    cm_dict = {}

    for predictor in predictors:
        if predictor not in data.columns:
            raise KeyError(f"Predictor {predictor} not found in the data.")
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
        
        cm_dict[predictor] = [int(tn), int(fp), int(fn), int(tp)]

    return cm_dict

###############

def main():
    # argument parsing
    parser = argparse.ArgumentParser(description="Evaluate predictors using PR AUC and/or ROC AUC")
    parser.add_argument('--ifile', help="Input file containing the prediction scores in TSV format (default: STDIN)", default=None)
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

    # get confusion matrix
    cm = get_confusion_matrix(data, args.predictors)

    # write the results to the output fileconfusion_matrix
    with open(args.ofile, 'w') as ofile:
        ofile.write(f"Tool\tTN\tFP\tFN\tTP\n")
        for predictor in args.predictors:
            ofile.write(f"{predictor}\t{cm[predictor][0]}\t{cm[predictor][1]}\t{cm[predictor][2]}\t{cm[predictor][3]}\n")

if __name__ == "__main__":
    main()
