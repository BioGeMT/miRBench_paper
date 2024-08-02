import argparse
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import numpy as np
import os

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_seed_match(target, mirna):
    rc_target = reverse_complement(target)
    seed6, seed7, seed8 = mirna[1:7], mirna[1:8], mirna[1:9]
    if seed8 in rc_target: return '8mer'
    elif seed7 in rc_target: return '7mer'
    elif seed6 in rc_target: return '6mer'
    return 'None'

def process_dataset(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df['seed_type'] = df.apply(lambda row: find_seed_match(row['seq.g'], row['seq.m']), axis=1)
    seed_percentages = df['seed_type'].value_counts(normalize=True) * 100
    return seed_percentages[['6mer', '7mer', '8mer']]


def plot_pooled_seed_prevalence(datasets, output_file):
    fig, ax = plt.subplots(figsize=(7, 5))  # Even more compact figure size
    
    data = []
    labels = []
    for file_path in datasets:
        seed_percentages = process_dataset(file_path)
        data.append(seed_percentages.values)
        labels.append(os.path.splitext(os.path.basename(file_path))[0])
    
    bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.5)  # Reduced box width further
    
    color = '#FFA500'  # Orange
    for patch in bp['boxes']:
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_linewidth(1)
    
    # Customize whiskers and caps
    for whisker in bp['whiskers']:
        whisker.set(color='black', linewidth=1)
    for cap in bp['caps']:
        cap.set(color='black', linewidth=1)
    
    # Customize medians
    for median in bp['medians']:
        median.set(color='black', linewidth=1.5)
    
    ax.set_xticklabels(labels, fontsize=12, fontweight='bold', rotation=0)
    ax.set_ylabel('% seed', fontsize=14, fontweight='bold')
    ax.set_title('Pooled Seed Prevalence Distribution', fontsize=16, fontweight='bold')
    
    ax.tick_params(axis='both', which='major', labelsize=12, labelcolor='black')
    
    ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.5)
    ax.set_axisbelow(True)
    
    # Set y-axis limits and ticks
    ax.set_ylim(0, 18)
    ax.set_yticks(range(0, 19, 2))
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Make left and bottom spines bold
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Adjust layout to reduce space between boxes
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate pooled seed prevalence distribution plot for multiple datasets.')
    parser.add_argument('datasets', nargs='+', help='Paths to datasets')
    parser.add_argument('-o', '--output', default='pooled_seed_prevalence_plot.png', help='Output file name (default: pooled_seed_prevalence_plot.png)')
    
    args = parser.parse_args()
    
    plot_pooled_seed_prevalence(args.datasets, args.output)

if __name__ == "__main__":
    main()