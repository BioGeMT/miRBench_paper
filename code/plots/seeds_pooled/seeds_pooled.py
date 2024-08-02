import argparse
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_seed_match(target, mirna):
    rc_target = reverse_complement(target)
    seed6, seed7, seed8 = mirna[1:7], mirna[1:8], mirna[1:9]
    if seed8 in rc_target: return '8mer'
    elif seed7 in rc_target: return '7mer'
    elif seed6 in rc_target: return '6mer'
    return 'None'

def process_datasets(datasets):
    all_data = []
    for name, file_path in datasets.items():
        df = pd.read_csv(file_path, sep='\t')
        df['seed_type'] = df.apply(lambda row: find_seed_match(row['seq.g'], row['seq.m']), axis=1)
        all_data.append(df)
    
    combined_df = pd.concat(all_data)
    mirna_counts = combined_df['miRNA_fam'].value_counts()
    top_mirnas = mirna_counts.nlargest(10).index.tolist()
    
    df_top = combined_df[combined_df['miRNA_fam'].isin(top_mirnas)]
    seed_percentages = df_top.groupby('miRNA_fam')['seed_type'].value_counts(normalize=True).unstack() * 100
    seed_percentages = seed_percentages.reindex(columns=['6mer', '7mer', '8mer']).fillna(0)
    
    return seed_percentages

def plot_seed_prevalence(seed_percentages, output_file):
    fig, ax = plt.subplots(figsize=(20, 12))  # Increased figure size further
    
    colors = ['#FFA500', '#90EE90', '#87CEFA']
    bar_width = 0.25
    index = range(len(seed_percentages.index))
    
    for i, seed_type in enumerate(['6mer', '7mer', '8mer']):
        values = seed_percentages[seed_type]
        ax.bar([x + i*bar_width for x in index], values, bar_width, color=colors[i], label=seed_type)
    
    ax.set_xlabel('miRNA Family', fontweight='bold', fontsize=20)  # Increased font size
    ax.set_ylabel('Percentage', fontweight='bold', fontsize=20)  # Increased font size
    ax.set_title('Seed Prevalence for Top 10 Expressed miRNAs', fontweight='bold', fontsize=24)  # Increased font size
    ax.set_xticks([x + bar_width for x in index])
    ax.set_xticklabels(seed_percentages.index, rotation=45, ha='right', fontsize=20, fontweight='bold')  # Increased font size to 20
    
    ax.set_ylim(0, 40)
    ax.set_yticks(range(0, 41, 5))
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    
    # Set tick parameters
    ax.tick_params(axis='both', which='major', labelsize=18)  # Increased y-axis label size
    
    # Manually set font weight for tick labels
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    legend = ax.legend(title='Seed Type', title_fontsize=20, fontsize=18, loc='upper right')  # Increased font sizes
    plt.setp(legend.get_title(), fontweight='bold')
    plt.setp(legend.get_texts(), fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate miRNA seed prevalence plot for three datasets.')
    parser.add_argument('dataset1', help='Path to first dataset')
    parser.add_argument('dataset2', help='Path to second dataset')
    parser.add_argument('dataset3', help='Path to third dataset')
    parser.add_argument('-o', '--output', default='seed_prevalence_plot.png', help='Output file name (default: seed_prevalence_plot.png)')
    
    args = parser.parse_args()
    
    datasets = {
        'Dataset1': args.dataset1,
        'Dataset2': args.dataset2,
        'Dataset3': args.dataset3
    }
    
    seed_percentages = process_datasets(datasets)
    plot_seed_prevalence(seed_percentages, args.output)

if __name__ == "__main__":
    main()