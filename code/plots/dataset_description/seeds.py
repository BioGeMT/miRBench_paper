import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio.Seq import Seq
import argparse
from collections import OrderedDict
import os
from seed_utils import find_seed_match


def process_dataset(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df['seed_type'] = df.apply(lambda row: find_seed_match(row['gene'], row['noncodingRNA']), axis=1)
    seed_counts = df['seed_type'].value_counts()
    return seed_counts

def create_tsv(datasets, tsv_output):
    data = []
    for name, file_path in datasets.items():
        seed_counts = process_dataset(file_path)
        data.append(pd.DataFrame({'Dataset': name, 'Seed': seed_counts.index, 'Count': seed_counts.values}))

    df = pd.concat(data)
    df_pivot = df.pivot(index='Dataset', columns='Seed', values='Count')
    df_pivot = df_pivot.reindex(columns=['Seed6mer', 'Seed7mer', 'Seed8mer', 'SeedNonCanonical', 'None'])
    df_pivot = df_pivot.fillna(0)

    # Reorder rows to match input order
    df_pivot = df_pivot.reindex(datasets.keys())

    # Calculate percentages
    df_pivot_pct = df_pivot.div(df_pivot.sum(axis=1), axis=0) * 100

    # Combine counts and percentages for TSV output
    combined_df = df_pivot.astype(int).join(df_pivot_pct, lsuffix='_count', rsuffix='_pct')

    # Save combined counts and percentages to TSV file
    combined_df.to_csv(tsv_output, sep='\t')
    print(f"TSV file saved as {tsv_output}")

    return combined_df

def plot_seed_prevalence(tsv_file, output_file):
    # Read the TSV file
    combined_df = pd.read_csv(tsv_file, sep='\t', index_col=0)

    # Separate count and percentage data
    df_pivot = combined_df.filter(regex='_count$').rename(columns=lambda x: x.replace('_count', ''))
    df_pivot_pct = combined_df.filter(regex='_pct$').rename(columns=lambda x: x.replace('_pct', ''))

    # Set figure size
    fig, ax = plt.subplots(figsize=(3500 / 300, 3500 / 300), dpi=300)

    # Set up the color palette (order as specified)
    colors = ['#662506', '#993404', '#cc4c02', '#ec7014', '#fe9929']  # Colorblind-safe orange/brown

    # Create the vertical stacked bar plot
    bar_width = 3.9
    bar_positions = range(3, 3 + 6 * len(df_pivot_pct), 6)

    bottom = [0] * len(df_pivot_pct)

    for i, seed in enumerate(['None', 'SeedNonCanonical', 'Seed6mer', 'Seed7mer', 'Seed8mer']):
        values = df_pivot_pct[seed].values
        ax.bar(bar_positions, values, bar_width, bottom=bottom, color=colors[i], edgecolor='none', linewidth=0.5)
        bottom = [b + v for b, v in zip(bottom, values)]

    # Add black borders to the entire columns
    for pos in bar_positions:
        ax.bar(pos, 100, bar_width, fill=False, edgecolor='none', linewidth=2.0)

    # Adjust y-axis
    ax.set_yticks(range(0, 101, 10))
    ax.set_yticklabels(["0%", "", "20%", "", "40%", "", "60%", "", "80%", "", "100%"],fontsize=30)
    #ax.set_yticklabels([f'{i}%' for i in range(0, 101, 10)], fontsize=30)
    ax.yaxis.grid(True, linestyle=':', alpha=0.8, color='black')

    # Increase width of x and y axis lines
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Adjust x-axis with reduced font size and specified order
    x_labels = df_pivot_pct.index
    plt.xticks(bar_positions, x_labels, fontsize=26, rotation=0, ha='center', weight='medium')

    # Move x-axis ticks down
    ax.xaxis.set_tick_params(pad=15)

    # Create custom square legend handles
    square_size = 24
    legend_handles = [Patch(facecolor=color, edgecolor='none') for color in colors]

    # Add legend to the bottom of the plot with corrected labels
    legend_labels = ['None', 'SeedNonCanonical', 'Seed6mer', 'Seed7mer', 'Seed8mer']
    legend = ax.legend(legend_handles, legend_labels, 
                       loc='upper center', bbox_to_anchor=(0.5, -0.10),
                       ncol=5, columnspacing=0.6, handletextpad=0.2, handlelength=1, handleheight=1,
                       prop={'size': 20, 'weight': 'medium'})

    # Adjust legend marker size and alignment
    for handle in legend.get_patches():
        handle.set_height(square_size)
        handle.set_width(square_size)

    # Center legend text with handles and reduce padding
    for text, handle in zip(legend.get_texts(), legend.get_patches()):
        text.set_ha('center')
        text.set_position((handle.get_x() + handle.get_width()/2, text.get_position()[1]))

    # Remove legend border
    legend.get_frame().set_edgecolor('none')

    # Adjust x-axis limits to extend beyond the legend width and add padding
    ax.set_xlim(0, len(df_pivot_pct) * 6)

    plt.tight_layout()

    # Save the plot with extra space at the bottom for the legend
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"Plot saved as {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate miRNA seed prevalence plot and TSV for datasets in specified order.')
    parser.add_argument('datasets', nargs='+', help='Paths to dataset files in the order you want them processed')
    parser.add_argument('-o', '--output', default='seed_prevalence_plot.png', help='Output file name for plot (default: seed_prevalence_plot.png)')
    parser.add_argument('-t', '--tsv', default='seed_prevalence_results.tsv', help='Output file name for TSV results (default: seed_prevalence_results.tsv)')

    args = parser.parse_args()

    # Retain dataset names, excluding file extensions
    datasets = OrderedDict()
    for file_path in args.datasets:
        dataset_name = os.path.splitext(os.path.basename(file_path))[0]
        datasets[dataset_name] = file_path

    # Create TSV file
    create_tsv(datasets, args.tsv)

    # Create plot using the TSV data
    plot_seed_prevalence(args.tsv, args.output)

if __name__ == "__main__":
    main()
