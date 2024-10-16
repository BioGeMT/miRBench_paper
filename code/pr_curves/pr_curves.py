import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import precision_recall_fscore_support
import argparse
import sys

def load_data(input_file):
    # load the data from the input file
    return pd.read_csv(input_file, sep='\t')

def plot_pr_curve(data, predictors, figsize=(6, 6), dpi=300, title=None):

    # set color scheme with 20 colors
    colors = plt.cm.tab20(np.linspace(0, 1, 11))
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)

    markers = ['o', 's', '^', 'X']
    i = 0

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    for predictor in predictors:

        if predictor not in data.columns:
            raise KeyError(f"Predictor {predictor} not found in the data.")
            
        if predictor.startswith('Seed'):
            p, r, _, _ = precision_recall_fscore_support(data['label'].values, data[predictor].values, average='binary')
            ax.plot(r, p, color='black', marker=markers[i], label=predictor)
            i += 1
        else:
            precision, recall, _ = precision_recall_curve(data['label'], data[predictor])
            ax.plot(recall, precision, label=predictor)
    
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    ax.legend()
    ax.set_title(title)
    fig.tight_layout()

    return fig, ax

def main():
    # argument parsing
    parser = argparse.ArgumentParser(description="Precision-Recall plots creator. ")
    parser.add_argument('--ifile', help="Input file (default: STDIN)", default=None)
    parser.add_argument('--predictors', help="List of predictor names (default: all)", default=None)
    parser.add_argument('--ofile', help="Output file (default: STDOUT)", default=None)
    parser.add_argument('--title', help="Title of the plot", default="Precision-Recall Curve")
    args = parser.parse_args()

    # if ifile is none, set it to sys.stdin
    if args.ifile is None:
        args.ifile = sys.stdin

    # load the data
    data = load_data(args.ifile)

    # if ofile is none, set it to sys.stdout
    if args.ofile is None:
        args.ofile = sys.stdout

    # if predictors is none, set it to all columns, except columns noncodingRNA, gene and label
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
            'Seed6merBulgeOrMismatch'] # in chronological order

    # plot precision-recall curves for all predictors
    fig, ax = plot_pr_curve(data, args.predictors, title=args.title)

    plt.savefig(args.ofile, format='svg')
    plt.close()

if __name__ == "__main__":
    main()
