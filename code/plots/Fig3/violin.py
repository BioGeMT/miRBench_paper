import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from collections import OrderedDict
import os
from seed_utils import find_seed_match

def process_dataset(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df['seed_type'] = df.apply(lambda row: find_seed_match(row['seq.g'], row['seq.m']), axis=1)
    seed_counts = df['seed_type'].value_counts()
    total = seed_counts.sum()
    seed_percentages = (seed_counts / total) * 100
    return seed_percentages

def create_combined_df(datasets):
    data = []
    for name, file_path in datasets.items():
        seed_percentages = process_dataset(file_path)
        data.append(pd.DataFrame({'Dataset': name, 'Seed': seed_percentages.index, 'Percentage': seed_percentages.values}))

    df = pd.concat(data)
    df_pivot = df.pivot(index='Dataset', columns='Seed', values='Percentage')
    df_pivot = df_pivot.reindex(columns=['SeedNonCanonical', 'Seed6mer', 'Seed7mer', 'Seed8mer'])
    df_pivot = df_pivot.fillna(0)

    # reorder rows to match input order
    df_pivot = df_pivot.reindex(datasets.keys())

    return df_pivot

def plot_seed_prevalence_violin(df_pivot, output_file):
    # set up the plot with adjusted dimensions
    fig, ax = plt.subplots(figsize=(3500 / 300, 3500 / 300), dpi=300)

    # define color palette with the specified colors
    colors = ['#f0e442', '#0072b2', '#cc79a7']

    # create the violin plot with borders
    positions = np.arange(len(df_pivot.index))
    violins = []

    for i, dataset in enumerate(df_pivot.index):
        data = df_pivot.loc[dataset]
        parts = ax.violinplot(data, positions=[i], showmeans=False, showmedians=False, showextrema=False,
                              widths=0.7)  
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i % len(colors)])
            pc.set_edgecolor('black')
            pc.set_linewidth(1.5)  
            pc.set_alpha(1)
        violins.append(parts['bodies'][0])

    # add lines for 25th, 50th, and 75th percentiles
    for i, violin in enumerate(violins):
        data = df_pivot.iloc[i]
        percentiles = np.percentile(data, [25, 50, 75])

        # get the path of the violin plot
        path = violin.get_paths()[0]
        vertices = path.vertices

        for j, percentile in enumerate(percentiles):
            # find the width of the violin at this percentile
            idx = np.argmin(np.abs(vertices[:, 1] - percentile))
            x_left = vertices[idx, 0]
            x_right = vertices[-idx, 0]

            linewidth = 2.5 if j == 1 else 2  # bolder lines for all, even bolder for median
            ax.hlines(percentile, x_left, x_right, color='black', linewidth=linewidth)

    # adjust y-axis to max 60% with 10% increments
    ax.set_ylim(0, 40)
    ax.set_yticks(range(0, 51, 10))
    ax.set_yticklabels([f'{i}%' for i in range(0, 51, 10)], fontsize=28)  # increased font size
    ax.yaxis.grid(True, linestyle=':', alpha=0.8, color='black')

    # increase width of x and y axis lines
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    # remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # adjust x-axis with horizontal ticks
    x_labels = df_pivot.index
    ax.set_xticks(positions)
    ax.set_xticklabels(x_labels, fontsize=24, rotation=0, ha='center', weight='semibold')  

    # move x-axis ticks down and add more padding
    ax.xaxis.set_tick_params(pad=15)

    # remove x and y labels
    ax.set_xlabel('')
    ax.set_ylabel('')

    # adjust the plot limits to reduce space around violins
    ax.set_xlim(-0.5, len(positions) - 0.5)

    # add padding around the plot
    plt.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.15)  

    plt.tight_layout()

    # save the plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight', format='png')
    print(f"plot saved as {output_file}")

def parse_arguments():
    parser = argparse.ArgumentParser(description='generate mirna seed prevalence violin plot for datasets in specified order.')
    parser.add_argument('datasets', nargs='+', help='paths to dataset files in the order you want them processed')
    parser.add_argument('-o', '--output', default='seed_prevalence_violin_plot.png', help='output file name for plot (default: seed_prevalence_violin_plot.png)')
    parser.add_argument('-t', '--tsv', default='seed_prevalence_data.tsv', help='output file name for tsv data (default: seed_prevalence_data.tsv)')
    return parser.parse_args()

def create_datasets_dict(dataset_paths):
    datasets = OrderedDict()
    for file_path in dataset_paths:
        dataset_name = os.path.splitext(os.path.basename(file_path))[0]
        datasets[dataset_name] = file_path
    return datasets

def save_data_to_tsv(df, file_path):
    df.to_csv(file_path, sep='\t')
    print(f"data saved as {file_path}")

def main():
    args = parse_arguments()
    datasets = create_datasets_dict(args.datasets)
    df_pivot = create_combined_df(datasets)
    save_data_to_tsv(df_pivot, args.tsv)
    plot_seed_prevalence_violin(df_pivot, args.output)

if __name__ == "__main__":
    main()