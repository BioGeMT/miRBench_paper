import json
import argparse
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve


COLOR_PALETTE = ['#E69F00','#56B4E9','#009E73','#F0E442','#CC79A7','#0072B2','#D55E00', '#000000']


def generate_random_predictions(length):
    return [random.uniform(0, 1) for _ in range(length)]


def calculate_scanning_length_at_sensitivity(predictions, labels, desired_sensitivity):
    """
    Function to calculate average scanning length (1/FPR) at a given sensitivity threshold.
    It calculates the FPR corresponding to the desired sensitivity and returns 1/FPR.
    """
    fpr, tpr, thresholds = roc_curve(labels, predictions)
    
    # Find the threshold corresponding to the desired sensitivity
    # Since TPR is sorted in increasing order, we can interpolate or find the closest one
    indices = np.where(tpr >= desired_sensitivity)[0]
    if len(indices) == 0:
        # Desired sensitivity not achievable
        print(f"Desired sensitivity {desired_sensitivity} not achievable.")
        return np.nan  # Or 0 or some other value

    index = indices[0]  # First index where TPR >= desired_sensitivity
    threshold = thresholds[index]
    fpr_at_threshold = fpr[index]

    if fpr_at_threshold == 0:
        scanning_length = np.inf
    else:
        scanning_length = 1 / fpr_at_threshold  # Scanning length is the inverse of FPR
    
    return scanning_length


def plot_scanning_length_based_on_sensitivity(df, labels_column, methods, sensitivities, max_value=2000, dpi=300, font_size=9, title='', show_y_axis_text=True, output_json_path='output/scanning_length_results.json'):
    """
    This function computes the scanning length (1/FPR) for each method at desired sensitivities and plots it.
    
    Parameters:
    - df: DataFrame containing the predictions and labels
    - labels_column: The name of the column containing the true labels
    - methods: A list of column names for the methods to be plotted
    - sensitivities: A list of desired sensitivities for each method
    - max_value: Maximum value to cap the infinite values, default is 2000 bp
    """
    labels = df[labels_column]
    
    # Calculate scanning lengths for each method
    lengths = []
    result={}
    for method, sensitivity in zip(methods, sensitivities):
        predictions = df[method]
        length = calculate_scanning_length_at_sensitivity(predictions, labels, sensitivity)
        if np.isinf(length) or np.isnan(length):
            length = max_value  # Replace infinite or NaN values with a large finite number
        lengths.append(length)
        result[f"{method}_{sensitivity}"]=length
        
    with open(output_json_path, 'w') as f:
        json.dump(result, f, indent=4)

    # Create the figure and axis for the plot
    if show_y_axis_text:
        fig, ax = plt.subplots(figsize=(8, 4), dpi=dpi)
    else:
        fig, ax = plt.subplots(figsize=(5, 3), dpi=dpi)

    # Bar positions and labels
    y_pos = np.arange(len(methods))

    # Add small space between y-axis and bars
    small_offset = lengths[-1] * 0.01  # X% of the last bar length as offset

    # Plotting the bars with a left offset
    ax.barh(
        y_pos, lengths, align='center', left=small_offset, 
        color=COLOR_PALETTE
    )

    # Adjust x-axis limits to accommodate the offset
    ax.set_xlim(0, max(lengths) + small_offset * 2)

    # Add text labels
    for i, v in enumerate(lengths):
        ax.text(v + small_offset, i, f"{v:.0f} bp", va='center', fontsize=font_size)

    ax.set_xticks([])
    if show_y_axis_text:
        ax.set_yticks(y_pos)
    else:
        ax.set_yticks([])
    ax.set_yticklabels([f"{meth} (Sens={round(sens, 3)})" for meth, sens in zip(methods, sensitivities)])
    ax.set_xlabel('Average Scanning Length (1/FPR) [bp] at recall 0.5')
    ax.set_title(title)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.tight_layout()


def main(args):
    dataset_path = args.dataset_dir + args.dataset_name
    
    predictions = pd.read_csv(f'{dataset_path}.tsv', sep='\t', header=0)
    predictions = predictions.rename(args.manakov_rename_map, axis=1)
    predictions["random"] = generate_random_predictions(predictions.shape[0])
    
    plot_scanning_length_based_on_sensitivity(
        predictions, labels_column='label', methods=args.method_names, 
        sensitivities=args.thresholds,
        show_y_axis_text=False,
        output_json_path='output/scanning_length_results.json',
    )
    
    title_suffix = ".avrg_fp_per_kbp_sensitivity_threshold"
    plt.savefig(f"output/{args.dataset_name}{title_suffix}.svg", format='svg')
    plt.savefig(f"output/{args.dataset_name}{title_suffix}.png", format='png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Average Scanning Length per Sensitivity Threshold Plot')
    parser.add_argument('--dataset-dir', type=str, default='predictions/', 
                       help='Directory containing the dataset')
    parser.add_argument('--output-dir', type=str, default='output/', 
                       help='Directory containing the figures and intermediate results')
    parser.add_argument('--dataset-name', type=str, 
                       default='AGO2_eCLIP_Manakov2022_100_CNN_predictions', 
                       help='Name of the dataset file')
    parser.add_argument('--thresholds', type=float, nargs='+', 
                       default=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 
                       help='Thresholds for scanning length plots')
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