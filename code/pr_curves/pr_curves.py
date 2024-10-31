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

def plot_pr_curve(data, predictors, title, dpi):

    # set color scheme with 20 colors
    #colors = plt.cm.tab20(np.linspace(0, 1, 11))
    colors = ['#cc79a7', '#d55e00', '#0072b2', '#f0e442', '#009e73', '#56b4e9', '#e69f00'] # 000000
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)

    markers = ['o', 's', '^', 'X']
    i = 0

    fig, ax = plt.subplots(figsize=(6, 6), dpi=int(dpi))

    for predictor in predictors:

        if predictor not in data.columns:
            raise KeyError(f"Predictor {predictor} not found in the data.")
            
        if predictor.startswith('Seed'):
            p, r, _, _ = precision_recall_fscore_support(data['label'].values, data[predictor].values, average='binary')
            ax.plot(r, p, color='black', marker=markers[i], markersize=12, label=predictor)
            i += 1
        else:
            precision, recall, _ = precision_recall_curve(data['label'], data[predictor])
            ax.plot(recall, precision, label=predictor, linewidth=2.5)
    
        ax.set_aspect(aspect=1)

        ax.set_xlabel('Recall', fontsize=17, labelpad=-15)
        ax.set_ylabel('Precision', fontsize=18, labelpad=-10)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

        ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.tick_params(axis='both', length=8)

        ax.set_xticklabels(['0', '', '', '', '', '1'], fontsize=16)
        ax.set_yticklabels(['0', '', '', '', '', '1'], fontsize=16)

        ax.text(-0.05, -0.05, '0', fontsize=16, ha='center', va='center', transform=ax.transAxes)


    ax.legend()

    #Title at the top
    ax.set_title(title, fontsize=22, fontweight='bold', pad=10)
    # Title on the left, parallel to the y-axis
    ax.text(-0.15, 0.5, 'dataset_name', va='center', ha='center', rotation=90, fontsize=22, fontweight='bold', transform=ax.transAxes)
    # Adjust plot to fit the side text
    plt.subplots_adjust(left=0.2)

    fig.tight_layout()

    return fig, ax

def main():
    # argument parsing
    parser = argparse.ArgumentParser(description="Precision-Recall plots creator. ")
    parser.add_argument('--ifile', help="Input file (default: STDIN)", default=None)
    parser.add_argument('--predictors', help="List of predictor names (default: all)", default=None)
    parser.add_argument('--ofile', help="Output file (default: STDOUT)", default=None)
    parser.add_argument('--title', help="Title of the plot", default=None)
    parser.add_argument('--dpi', help="DPI (default: 300)", default=300)

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
    fig, ax = plot_pr_curve(data, args.predictors, title=args.title, dpi=args.dpi)

    plt.savefig(args.ofile)
    plt.close()

if __name__ == "__main__":
    main()
