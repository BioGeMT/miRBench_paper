import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from avrg_fp_per_sensitivity import COLOR_PALETTE, generate_random_predictions
from sklearn.metrics import precision_recall_curve


def plot_precision_at_recall(
    df, model_names, sensitivity_thresholds, labels, dpi=300, font_size=9, legend_font_size=9, title=''
):
    """
    Plots precision at specified recall (sensitivity) thresholds for the given models.

    Parameters:
    - df: DataFrame containing the predictions and true labels.
          It should have columns for model predictions and one column for true labels.
    - model_names: List of model prediction column names to plot.
    - sensitivity_thresholds: List of recall (sensitivity) thresholds to use.
    - labels: Column name for the true labels.
    """
    true_labels = df[labels]
    
    precision_at_sensitivities = []

    for model in model_names:
        precisions = []
        predictions = df[model]
        
        precision, recall, thresholds = precision_recall_curve(true_labels, predictions)
        
        for threshold in sensitivity_thresholds:
            idx_closest = np.argmin(np.abs(recall - threshold))
            precisions.append(precision[idx_closest])
        
        precision_at_sensitivities.append(precisions)
    
    precision_values = np.array(precision_at_sensitivities)

    n_models = len(model_names)
    n_sensitivities = len(sensitivity_thresholds)
    total_bar_width = 1  # Total width of all bars at one sensitivity threshold
    bar_width = total_bar_width / n_models  # Width of each bar
    group_spacing = 1.1  # Adjust this value to bring groups closer together (less than 1)
    index = np.arange(n_sensitivities) * group_spacing

    offset = (group_spacing - total_bar_width) / 2

    fig, ax = plt.subplots(figsize=(14, 5), dpi=dpi)
    
    small_offset = 0.1
    
    for i, model in enumerate(model_names):
        bars = ax.bar(
            index + offset + i * bar_width,
            precision_values[i],
            bar_width,
            label=model,
            color=COLOR_PALETTE[i % len(COLOR_PALETTE)]
        )
        
        for bar in bars:
            height = bar.get_height()
            ax.annotate(
                f'{height:.2f}',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # Slightly above the bar
                textcoords="offset points",
                ha='center',
                va='bottom',
                fontsize=font_size,
            )
    
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title(title)
    ax.set_xticks(index + group_spacing / 2)
    ax.set_xticklabels([f'{t}' for t in sensitivity_thresholds])
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.subplots_adjust(bottom=0.2, left=0.35, right=0.65)

    ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.15),  # Position below the plot
        ncol=len(model_names),  # Place all items in one row
        prop={'size': legend_font_size},
        frameon=False,
    )
    
    ax.set_ylim(bottom=-0.01)
    ax.set_ylim(-0.01, 1)
    
    plt.tight_layout()
    

def main(args):
    dataset_path = args.dataset_dir + args.dataset_name
    
    predictions = pd.read_csv(f'{dataset_path}.tsv', sep='\t', header=0)
    predictions = predictions.rename(args.manakov_rename_map, axis=1)
    predictions["random"] = generate_random_predictions(predictions.shape[0])
    
    plot_precision_at_recall(
        predictions, args.method_names, args.sensitivity_thresholds, labels='label',
    )
    
    title_suffix = ".precision_per_sensitivity_threshold"
    plt.savefig(f"output/{args.dataset_name}{title_suffix}.svg", format='svg')
    plt.savefig(f"output/{args.dataset_name}{title_suffix}.png", format='png')

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Precision per Sensitivity Threshold Plot')
    parser.add_argument('--dataset-dir', type=str, default='predictions/', 
                       help='Directory containing the dataset')
    parser.add_argument('--output-dir', type=str, default='output/', 
                       help='Directory containing the figures and intermediate results')
    parser.add_argument('--dataset-name', type=str, 
                       default='AGO2_eCLIP_Manakov2022_100_CNN_predictions', 
                       help='Name of the dataset file')
    parser.add_argument('--sensitivity-thresholds', type=float, nargs='+', 
                       default=[0.5, 0.33], 
                       help='Sensitivity thresholds for precision-recall plots')
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