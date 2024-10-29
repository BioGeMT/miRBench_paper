import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import sys
from matplotlib.patches import Patch

def load_data(file_path):
    try:
        return pd.read_csv(file_path or sys.stdin, sep='\t')
    except Exception as e:
        print(f"Error reading input data: {e}", file=sys.stderr)
        sys.exit(1)

def categorize_features(data):
    categories = ['five_prime_utr', 'three_prime_utr', 'intron', 'exon', 'Unknown']
    counts = {cat: 0 for cat in categories}
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

    return {cat: (count / total_entries) * 100 for cat, count in counts.items()}

def save_to_tsv(percentages, labels, output_file):
    pd.DataFrame(percentages, index=labels).to_csv(output_file, sep='\t', index_label='Dataset')

def create_plot(percentages, labels):
    categories = ["3' UTR", "5' UTR", 'Intron', 'Exon', 'Unknown']
    datasets = [
        [p['three_prime_utr'], p['five_prime_utr'], p['intron'], p['exon'], p['Unknown']]
        for p in percentages
    ]

    # create figure and axis
    fig, ax = plt.subplots(figsize=(8.27, 8.27))  # square figure, a4 width
    colors = ['#d55e00', '#e69f00', '#f0e442', '#0072b2', '#cc79a7']

    # set bar properties
    bar_width = 0.5
    bar_positions = [i * 0.75 for i in range(len(datasets))]

    # create stacked bars for each dataset
    for i, dataset in enumerate(datasets):
        bottom = 0
        for j, value in enumerate(dataset):
            ax.bar(bar_positions[i], value, bar_width, bottom=bottom, color=colors[j], edgecolor='black', linewidth=0.5)
            bottom += value

    # add black borders to the entire columns
    for pos in bar_positions:
        ax.bar(pos, 100, bar_width, fill=False, edgecolor='black', linewidth=1.5)

    # adjust y-axis
    ax.set_yticks(range(0, 101, 10))
    ax.set_yticklabels(["0%", "", "20%", "", "40%", "", "60%", "", "80%", "", "100%"],fontsize=19)
    ax.yaxis.grid(True, linestyle=':', alpha=0.8, color='black')

    # increase width of x and y axis lines
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    # remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # adjust x-axis with half-bold font
    x_labels = [label.replace('.tsv', '') for label in labels]
    plt.xticks(bar_positions, x_labels, fontsize=19, rotation=0, ha='center', weight='medium')

    # move x-axis ticks down by 1/3
    ax.xaxis.set_tick_params(pad=15)

    # create custom square legend handles
    square_size = 25  # reduced by 1/3 from 30
    legend_handles = [Patch(facecolor=color, edgecolor='black') for color in colors]

    # add legend below x-axis
    legend = ax.legend(legend_handles, categories, loc='upper center', bbox_to_anchor=(0.5, -0.10),
                       ncol=5, columnspacing=1, handlelength=1.5, handleheight=1.5, handletextpad=0.0,
                       prop={'size': 17, 'weight': 'medium'})

    # adjust legend marker size and alignment
    for handle in legend.get_patches():
        handle.set_height(square_size)
        handle.set_width(square_size)

    # center legend text with handles
    for text, handle in zip(legend.get_texts(), legend.get_patches()):
        text.set_ha('center')  # center text horizontally
        text.set_position((handle.get_x() + handle.get_width()/2, text.get_position()[1]))

    # remove legend border
    legend.get_frame().set_edgecolor('none')

    plt.tight_layout()
    return fig

def parse_arguments():
    parser = argparse.ArgumentParser(description="miRNA Feature Distribution Plot")
    parser.add_argument('datasets', nargs='*', help="Input dataset files")
    parser.add_argument('-o', '--output', help="Output file for plot (e.g., plot.png)")
    parser.add_argument('-t', '--tsv', help="Output file for TSV data (e.g., output.tsv)")
    return parser.parse_args()

def process_datasets(datasets):
    return [(load_data(file), os.path.basename(file)) for file in datasets]

def calculate_results(datasets):
    return [categorize_features(dataset) for dataset, _ in datasets]

def save_plot(fig, output_file):
    if output_file:
        fig.savefig(output_file, bbox_inches='tight', dpi=300)
    else:
        plt.show()

def main():
    args = parse_arguments()
    datasets = process_datasets(args.datasets)
    percentages = calculate_results(datasets)
    labels = [label for _, label in datasets]
    
    fig = create_plot(percentages, labels)
    save_plot(fig, args.output)
    if args.tsv:
        save_to_tsv(percentages, labels, args.tsv)

if __name__ == "__main__":
    main()