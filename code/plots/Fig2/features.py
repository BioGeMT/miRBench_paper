import pandas as pd
import matplotlib.pyplot as plt
import argparse
from collections import Counter
import os
import sys
from matplotlib.patches import Patch

def load_data(file_path):
    try:
        if file_path:
            return pd.read_csv(file_path, sep='\t')
        else:
            return pd.read_csv(sys.stdin, sep='\t')
    except Exception as e:
        print(f"Error reading input data: {e}", file=sys.stderr)
        sys.exit(1)

def categorize_features(data):
    categories = ['five_prime_utr', 'three_prime_utr', 'intron', 'exon', 'Unknown']
    counts = {cat: 0 for cat in categories}
    other_features = Counter()
    total_entries = len(data)
    
    for feature in data['feature']:
        feature = str(feature).lower()
        if 'five_prime_utr' in feature:
            counts['five_prime_utr'] += 1
        elif 'three_prime_utr' in feature:
            counts['three_prime_utr'] += 1
        elif 'intron' in feature:
            counts['intron'] += 1
        elif 'exon' in feature:
            counts['exon'] += 1
        else:
            counts['Unknown'] += 1
            other_features[feature] += 1
    
    percentages = {cat: (count / total_entries) * 100 for cat, count in counts.items()}
    return percentages, other_features

def create_plot(percentages1, percentages2, percentages3, label1, label2, label3):
    categories = ["5' UTR", "3' UTR", 'Intron', 'Exon', 'Unknown']
    dataset1_values = [percentages1['five_prime_utr'], percentages1['three_prime_utr'], percentages1['intron'], percentages1['exon'], percentages1['Unknown']]
    dataset2_values = [percentages2['five_prime_utr'], percentages2['three_prime_utr'], percentages2['intron'], percentages2['exon'], percentages2['Unknown']]
    dataset3_values = [percentages3['five_prime_utr'], percentages3['three_prime_utr'], percentages3['intron'], percentages3['exon'], percentages3['Unknown']]
    
    fig, ax = plt.subplots(figsize=(8.27, 8.27))  # Square figure, A4 width
    bottom1, bottom2, bottom3 = 0, 0, 0
    colors = ['#d55e00', '#e69f00', '#f0e442', '#0072b2', '#cc79a7']
    
    bar_width = 0.5
    bar_positions = [0, 0.75, 1.5]
    for i, cat in enumerate(categories):
        ax.bar(bar_positions[0], dataset1_values[i], bar_width, bottom=bottom1, color=colors[i], edgecolor='black', linewidth=0.5)
        ax.bar(bar_positions[1], dataset2_values[i], bar_width, bottom=bottom2, color=colors[i], edgecolor='black', linewidth=0.5)
        ax.bar(bar_positions[2], dataset3_values[i], bar_width, bottom=bottom3, color=colors[i], edgecolor='black', linewidth=0.5)
        bottom1 += dataset1_values[i]
        bottom2 += dataset2_values[i]
        bottom3 += dataset3_values[i]
    
    # Add black borders to the entire columns
    for pos in bar_positions:
        ax.bar(pos, 100, bar_width, fill=False, edgecolor='black', linewidth=1.5)
    
    # Adjust y-axis
    ax.set_yticks(range(0, 101, 10))
    ax.set_yticklabels([f'{i}%' for i in range(0, 101, 10)], fontsize=22)
    ax.yaxis.grid(True, linestyle=':', alpha=0.8, color='black')
    
    # Increase width of x and y axis lines
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Adjust x-axis with half-bold font
    x_labels = [label.replace('.tsv', '') for label in [label1, label2, label3]]
    plt.xticks(bar_positions, x_labels, fontsize=19, rotation=0, ha='center', weight='semibold')
    
    # Move x-axis ticks down by 1/3
    ax.xaxis.set_tick_params(pad=15)
    
    # Create custom square legend handles
    square_size = 20  # Reduced by 1/3 from 30
    legend_handles = [Patch(facecolor=color, edgecolor='black') for color in colors]
    
    # Add legend below x-axis
    legend = ax.legend(legend_handles, categories, loc='upper center', bbox_to_anchor=(0.45, -0.10),
                       ncol=5, columnspacing=1, handlelength=1.5, handleheight=1.5,
                       prop={'size': 16, 'weight': 'semibold'})
    
    # Adjust legend marker size and alignment
    for handle in legend.get_patches():
        handle.set_height(square_size)
        handle.set_width(square_size)
    
    # Center legend text with handles
    for text, handle in zip(legend.get_texts(), legend.get_patches()):
        text.set_ha('center')  # Center text horizontally
        text.set_position((handle.get_x() + handle.get_width()/2, text.get_position()[1]))
    
    # Remove legend border
    legend.get_frame().set_edgecolor('none')
    
    plt.tight_layout()
    return fig

def print_other_features(other_features, dataset_name):
    print(f"\nDetailed breakdown of 'Unknown' features for {dataset_name}:")
    total = sum(other_features.values())
    sorted_features = sorted(other_features.items(), key=lambda x: x[1], reverse=True)
    for feature, count in sorted_features:
        percentage = (count / total) * 100
        print(f"  {feature}: {count} ({percentage:.2f}%)")

def main():
    parser = argparse.ArgumentParser(description="miRNA Feature Distribution Plot")
    parser.add_argument('datasets', nargs='*', help="Input dataset files (up to 3)")
    parser.add_argument('-o', '--output', help="Output file for plot (e.g., plot.png)")
    args = parser.parse_args()

    if len(args.datasets) > 3:
        print("Error: Maximum of 3 datasets allowed.", file=sys.stderr)
        sys.exit(1)

    datasets = []
    labels = []

    # Load datasets
    for dataset_file in args.datasets:
        datasets.append(load_data(dataset_file))
        labels.append(os.path.basename(dataset_file))

    # Calculate percentages and get other features
    results = [categorize_features(dataset) for dataset in datasets]
    percentages = [r[0] for r in results]
    other_features = [r[1] for r in results]

    # Ensure we have exactly 3 datasets (pad with empty data if necessary)
    while len(percentages) < 3:
        percentages.append({cat: 0 for cat in ['five_prime_utr', 'three_prime_utr', 'intron', 'exon', 'Unknown']})
        labels.append("Empty")

    # Create plot
    fig = create_plot(percentages[0], percentages[1], percentages[2], labels[0], labels[1], labels[2])

    # Save plot if output file is specified
    if args.output:
        fig.savefig(args.output, bbox_inches='tight', dpi=300)
    else:
        plt.show()

    # Print percentages and other features for verification
    for label, perc, other in zip(labels, percentages, other_features):
        print(f"\n{label} percentages:")
        for cat, value in perc.items():
            print(f"  {cat}: {value:.2f}%")
        print_other_features(other, label)

if __name__ == "__main__":
    main()