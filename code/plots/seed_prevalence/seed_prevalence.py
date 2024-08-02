import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import argparse

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_seed_match(target, mirna):
    rc_target = reverse_complement(target)
    seed6 = mirna[1:7]
    seed7 = mirna[1:8]
    seed8 = mirna[1:9]
    
    if seed8 in rc_target:
        return '8mer'
    elif seed7 in rc_target:
        return '7mer'
    elif seed6 in rc_target:
        return '6mer'
    else:
        return 'none'

def process_dataset(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df['seed_type'] = df.apply(lambda row: find_seed_match(row['seq.g'], row['seq.m']), axis=1)
    seed_counts = df['seed_type'].value_counts(normalize=True) * 100
    return seed_counts

def plot_seed_prevalence(datasets, output_file):
    fig, ax = plt.subplots(figsize=(14, 8))  # Increased figure size
    
    data = []
    for name, file_path in datasets.items():
        seed_counts = process_dataset(file_path)
        data.append(pd.DataFrame({'Dataset': name, 'Seed': seed_counts.index, 'Percentage': seed_counts.values}))
    
    df = pd.concat(data)
    df_pivot = df.pivot(index='Dataset', columns='Seed', values='Percentage')
    df_pivot = df_pivot.reindex(columns=['none', '6mer', '7mer', '8mer'])
    df_pivot = df_pivot.fillna(0)
    
    # Set up the color palette
    colors = ['#FFA500', '#90EE90', '#87CEFA', '#FFC0CB']
    
    # Create the horizontal stacked bar plot
    df_pivot.plot(kind='barh', stacked=True, width=0.6, color=colors, ax=ax, legend=False)  # Reduced width to 0.6
    
    # Remove gridlines
    ax.grid(False)
    
    # Customize the plot
    ax.set_xlabel('Percentage', fontweight='bold', fontsize=18)  # Increased font size
    ax.set_ylabel('Dataset', fontweight='bold', fontsize=18)  # Increased font size
    
    # Bold tick labels with increased font size
    ax.tick_params(axis='both', which='major', labelsize=18)  # Increased font size
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Add a border to the plot
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.5)
    
    # Create a separate legend
    legend_labels = ['none', '6mer', '7mer', '8mer']
    legend_handles = [plt.Rectangle((0,0),1,1, facecolor=colors[i]) for i in range(len(legend_labels))]
    legend = fig.legend(legend_handles, legend_labels, 
                        loc='upper center', 
                        bbox_to_anchor=(0.5, 1.02),  # Moved closer to the plot
                        ncol=4, 
                        frameon=True, 
                        title='Seed Type',
                        title_fontsize=14,  # Increased font size
                        fontsize=12)  # Increased font size
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_linewidth(1)
    
    # Set the plot title
    fig.suptitle('Seed Prevalence per Dataset', fontweight='bold', fontsize=16, y=1.05)  # Increased font size and moved down
    
    # Adjust layout
    plt.tight_layout()
    
    # Adjust the subplot to leave less space at the top
    plt.subplots_adjust(top=0.9)  # Reduced top margin
    
    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate miRNA seed prevalence plot for three datasets.')
    parser.add_argument('Hejret2023', help='Path to Hejret2023 dataset')
    parser.add_argument('Klimentova2022', help='Path to Klimentova2022 dataset')
    parser.add_argument('Manakov2022', help='Path to Manakov2022 dataset')
    parser.add_argument('-o', '--output', default='seed_prevalence_plot.png', help='Output file name (default: seed_prevalence_plot.png)')
    
    args = parser.parse_args()
    
    datasets = {
        'Hejret2023': args.Hejret2023,
        'Klimentova2022': args.Klimentova2022,
        'Manakov2022': args.Manakov2022
    }
    
    plot_seed_prevalence(datasets, args.output)

if __name__ == "__main__":
    main()