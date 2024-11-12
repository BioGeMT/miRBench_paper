import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, auc
from plot_avrg_fp_per_sensitivity import COLOR_PALETTE
import argparse


def generate_random_predictions(length):
    return [random.uniform(0, 1) for _ in range(length)]


def calculate_auc_pr(predictions, method_names, label_column='label'):
    """Calculates AUC-PR for specified method columns in the predictions DataFrame."""
    auc_values = {}
    for method in method_names:
        if method in predictions.columns and pd.api.types.is_numeric_dtype(predictions[method]):
            precision, recall, _ = precision_recall_curve(predictions[label_column], predictions[method])
            auc_values[method] = auc(recall, precision)
        else:
            print(f"Warning: Method {method} not found or is not numeric.")
    return auc_values


def plot_auc_vs_dataset_size(auc_values, dataset_sizes, color_palette):
    """Plots AUC-PR against dataset size using a specified color palette."""
    x = [dataset_sizes[method] for method in auc_values.keys()]
    y = [auc_values[method] for method in auc_values.keys()]

    plt.figure(figsize=(5.5, 6))
    for i, method in enumerate(auc_values.keys()):
        plt.scatter(x[i], y[i], color=color_palette[i % len(color_palette)], label=method)
        plt.annotate(method, (x[i], y[i]))

    plt.xlabel("Dataset Size (log10)")
    plt.ylabel("AUC-PR")
    plt.annotate(method, (x[i] - 0.05, y[i]))
    plt.ylim(-0.05, 1)
    plt.grid(False)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    
def main(args):
    dataset_path = args.dataset_dir + args.dataset_name
    predictions = pd.read_csv(f'{dataset_path}.tsv', sep='\t', header=0)
    predictions = predictions.rename(args.manakov_rename_map, axis=1)
    predictions["random"] = generate_random_predictions(predictions.shape[0])
    
    # Verify label column exists
    if 'label' not in predictions.columns:
        raise ValueError("The dataset must contain a 'label' column for ground truth.")

    # Calculate AUC-PR for each specified method in method_names
    auc_values = calculate_auc_pr(predictions, args.method_names, label_column='label')
    
    # Dataset sizes (log10 scale for plotting)
    dataset_sizes = {
        'random': np.log10(predictions.shape[0]),
        'Manakov_subset_200': np.log10(200),
        'Manakov_subset_2k': np.log10(2000),
        'Manakov_subset_7720': np.log10(7720),
        'Manakov_subset_20k': np.log10(20000),
        'Manakov_subset_200k': np.log10(200000),
        'Manakov_2,524,246': np.log10(2524246)
    }
    
    plot_auc_vs_dataset_size(auc_values, dataset_sizes, COLOR_PALETTE)

    title_suffix = ".auc_pr_vs_train_data_size"
    plt.savefig(f"output/{args.dataset_name}{title_suffix}.svg", format='svg')
    plt.savefig(f"output/{args.dataset_name}{title_suffix}.png", format='png')

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='AUC-PR vs Dataset Size Plot')
    parser.add_argument('--dataset-dir', type=str, default='predictions/', 
                       help='Directory containing the dataset')
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
