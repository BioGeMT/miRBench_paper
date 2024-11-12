import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, auc
from sklearn.metrics import precision_recall_fscore_support
from plot_avrg_fp_per_sensitivity import COLOR_PALETTE, generate_random_predictions


def plot_pr_curve(
    data, predictors, figsize=(6, 6), dpi=300, font_size=12, title=None, output_json_path='output/pr_curve_results.json',
    colors = COLOR_PALETTE
):
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)

    markers = ['o', 's', '^', 'X']
    i = 0

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    result={}

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
            # print(f"{predictor} : {auc(recall, precision)}")
            result[predictor] = auc(recall, precision)
    
        ax.set_xlabel('Recall', fontsize=font_size)
        ax.set_ylabel('Precision', fontsize=font_size)

        # Set font size for tick labels
        ax.tick_params(axis='both', which='major', labelsize=font_size)
    # ax.legend()
    ax.set_title(title)
    fig.tight_layout()
    
    with open(output_json_path, 'w') as f:
        json.dump(result, f, indent=4)

    return fig, ax


def main(args):
    dataset_path = args.dataset_dir + args.dataset_name
    
    predictions = pd.read_csv(f'{dataset_path}.tsv', sep='\t', header=0)
    predictions = predictions.rename(args.manakov_rename_map, axis=1)
    predictions["random"] = generate_random_predictions(predictions.shape[0])
    
    fig, ax = plot_pr_curve(
        predictions,
        args.method_names,
        title="",
        output_json_path='output/pr_curve_results.json',
    )
    
    title_suffix = ".pr_curve"
    plt.savefig(f"output/{args.dataset_name}{title_suffix}.svg", format='svg')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Precision-Recall Plot')
    parser.add_argument('--dataset-dir', type=str, default='predictions/', 
                       help='Directory containing the dataset')
    parser.add_argument('--output-dir', type=str, default='output/', 
                       help='Directory containing the figures and intermediate results')
    parser.add_argument('--dataset-name', type=str, 
                       default='AGO2_eCLIP_Manakov2022_100_CNN_predictions', 
                       help='Name of the dataset file')
    parser.add_argument('--manakov-rename-map', type=lambda x: eval(x), default={
        'CNN_Manakov_full': 'Manakov_2,524,246',
        "CNN_Manakov_subset_200": 'Manakov_subset_200',
        "CNN_Manakov_subset_2k": 'Manakov_subset_2k',
        "CNN_Manakov_subset_7720": 'Manakov_subset_7720',
        "CNN_Manakov_subset_20k": 'Manakov_subset_20k',
        "CNN_Manakov_subset_200k": 'Manakov_subset_200k',
    }, help='Mapping for renaming Manakov dataset columns')
    parser.add_argument('--method-names', type=lambda x: eval(x), default=[
        "random", "Manakov_subset_200", "Manakov_subset_2k", "Manakov_subset_7720", 
        "Manakov_subset_20k", "Manakov_subset_200k", "Manakov_2,524,246"
    ], help='List of method names for plotting')
    args = parser.parse_args()
    main(args)