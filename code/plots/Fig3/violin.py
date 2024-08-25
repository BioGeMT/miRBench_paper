import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from Bio.Seq import Seq
from collections import OrderedDict
import os

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_seed_match(target, mirna):
    rc_target = reverse_complement(target)
    seed6 = mirna[1:7]
    seed7 = mirna[1:8]
    seed8 = 'A' + mirna[1:8]
    
    # Seed8mer
    if seed8 in rc_target:
        return 'Seed8mer'
    
    # Seed7mer
    if seed7 in rc_target or ('A' + seed6) in rc_target:
        return 'Seed7mer'
    
    # Seed6mer (including 6mer_m8/A1)
    if seed6 in rc_target or mirna[2:8] in rc_target or ('A' + mirna[1:6]) in rc_target:
        return 'Seed6mer'
    
    # SeedNonCanonical (formerly 6mer with bulge or mismatch)
    for pos in range(1, 7):
        for nt in ['A', 'C', 'G', 'T']:
            # bulges
            if (seed6[:pos] + nt + seed6[pos:]) in rc_target:
                return 'SeedNonCanonical'
            # mismatches
            if (seed6[:pos] + nt + seed6[pos+1:]) in rc_target:
                return 'SeedNonCanonical'
    
    # Check for 6mer_m8 and 6mer_A1 with bulge or mismatch
    seed6_m8 = mirna[2:8]
    seed6_A1 = 'A' + mirna[1:6]
    for seed in [seed6_m8, seed6_A1]:
        for pos in range(len(seed)):
            for nt in ['A', 'C', 'G', 'T']:
                # bulges
                if (seed[:pos] + nt + seed[pos:]) in rc_target:
                    return 'SeedNonCanonical'
                # mismatches
                if (seed[:pos] + nt + seed[pos+1:]) in rc_target:
                    return 'SeedNonCanonical'
    
    return 'none'

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
    df_pivot = df_pivot.reindex(columns=[ 'SeedNonCanonical', 'Seed6mer', 'Seed7mer', 'Seed8mer'])
    df_pivot = df_pivot.fillna(0)
    
    # Reorder rows to match input order
    df_pivot = df_pivot.reindex(datasets.keys())
    
    return df_pivot

def plot_seed_prevalence_violin(df_pivot, output_file):
    # Set up the plot with adjusted dimensions
    fig, ax = plt.subplots(figsize=(3500 / 300, 3500 / 300), dpi=300)
    
    # Define color palette with the specified colors
    colors = ['#f0e442', '#0072b2', '#cc79a7']
    
    # Create the violin plot with slightly narrower violins and thicker borders
    positions = np.arange(len(df_pivot.index))
    violins = []
    
    for i, dataset in enumerate(df_pivot.index):
        data = df_pivot.loc[dataset]
        parts = ax.violinplot(data, positions=[i], showmeans=False, showmedians=False, showextrema=False,
                              widths=0.7)  # Reduced width of violins
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i % len(colors)])
            pc.set_edgecolor('black')
            pc.set_linewidth(1.5)  # Set border width to 1.5
            pc.set_alpha(1)
        violins.append(parts['bodies'][0])
    
    # Add lines for 25th, 50th, and 75th percentiles
    for i, violin in enumerate(violins):
        data = df_pivot.iloc[i]
        percentiles = np.percentile(data, [25, 50, 75])
        
        # Get the path of the violin plot
        path = violin.get_paths()[0]
        vertices = path.vertices
        
        for j, percentile in enumerate(percentiles):
            # Find the width of the violin at this percentile
            idx = np.argmin(np.abs(vertices[:, 1] - percentile))
            x_left = vertices[idx, 0]
            x_right = vertices[-idx, 0]
            
            linewidth = 2.5 if j == 1 else 2  # Bolder lines for all, even bolder for median
            ax.hlines(percentile, x_left, x_right, color='black', linewidth=linewidth)
    
    # Adjust y-axis to max 60% with 10% increments
    ax.set_ylim(0, 40)
    ax.set_yticks(range(0, 51, 10))
    ax.set_yticklabels([f'{i}%' for i in range(0, 51, 10)], fontsize=28)  # Increased font size
    ax.yaxis.grid(True, linestyle=':', alpha=0.8, color='black')
    
    # Increase width of x and y axis lines
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Adjust x-axis with horizontal ticks
    x_labels = df_pivot.index
    ax.set_xticks(positions)
    ax.set_xticklabels(x_labels, fontsize=24, rotation=0, ha='center', weight='semibold')  # Increased font size
    
    # Move x-axis ticks down and add more padding
    ax.xaxis.set_tick_params(pad=15)
    
    # Remove x and y labels
    ax.set_xlabel('')
    ax.set_ylabel('')
    
    # Adjust the plot limits to reduce space around violins
    ax.set_xlim(-0.5, len(positions) - 0.5)
    
   # Add padding around the plot
    plt.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.15)  # Proper padding setup

    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight', format='png')
    print(f"Plot saved as {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate miRNA seed prevalence violin plot for datasets in specified order.')
    parser.add_argument('datasets', nargs='+', help='Paths to dataset files in the order you want them processed')
    parser.add_argument('-o', '--output', default='seed_prevalence_violin_plot.png', help='Output file name for plot (default: seed_prevalence_violin_plot.png)')
    
    args = parser.parse_args()
    
    # Retain dataset names, excluding file extensions
    datasets = OrderedDict()
    for file_path in args.datasets:
        dataset_name = os.path.splitext(os.path.basename(file_path))[0]
        datasets[dataset_name] = file_path

    # Process datasets and create combined dataframe
    df_pivot = create_combined_df(datasets)
    
    # Create violin plot
    plot_seed_prevalence_violin(df_pivot, args.output)

if __name__ == "__main__":
    main()
