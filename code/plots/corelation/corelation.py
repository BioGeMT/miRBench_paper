import sys
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

def load_and_process_data(file_path):
    df = pd.read_csv(file_path, sep='\t')
    mirna_counts = df['miRNA_fam'].value_counts()
    total_count = mirna_counts.sum()
    mirna_percentages = (mirna_counts / total_count) * 100
    return mirna_percentages

def calculate_correlations(datasets, file_names):
    all_mirnas = set()
    for dataset in datasets:
        all_mirnas.update(dataset.index)
    
    data = {os.path.splitext(name)[0]: dataset.reindex(all_mirnas).fillna(0) 
            for name, dataset in zip(file_names, datasets)}
    df = pd.DataFrame(data)
    
    corr_matrix = df.corr(method='pearson')
    return corr_matrix

def create_heatmap(corr_matrix, output_file):
    plt.figure(figsize=(12, 10))  # Increased figure size
    mask = np.tril(np.ones_like(corr_matrix, dtype=bool), k=-1)
    cmap = sns.color_palette("YlOrRd", as_cmap=True)
    
    # Create heatmap with larger font size for cell values
    ax = sns.heatmap(corr_matrix, mask=mask, annot=True, fmt=".2f", cmap=cmap, 
                     vmin=0, vmax=1, square=True, cbar_kws={'label': 'Correlation'},
                     annot_kws={'size': 14})  # Increased font size for cell values
    
    # Make labels bold and increase font size
    ax.set_xticklabels(ax.get_xticklabels(), weight='bold', size=12)
    ax.set_yticklabels(ax.get_yticklabels(), weight='bold', size=12)
    
    plt.title('miRNA Expression Correlation Heatmap', weight='bold', size=16)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('Correlation', size=14, weight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def main(input_files, output_tsv, output_heatmap):
    datasets = [load_and_process_data(file) for file in input_files]
    file_names = [os.path.basename(file) for file in input_files]
    corr_matrix = calculate_correlations(datasets, file_names)
    
    # Ensure diagonal is exactly 1
    np.fill_diagonal(corr_matrix.values, 1)
    
    corr_matrix.to_csv(output_tsv, sep='\t')
    create_heatmap(corr_matrix, output_heatmap)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py input1.tsv input2.tsv input3.tsv output_correlations.tsv output_heatmap.png")
        sys.exit(1)
    
    input_files = sys.argv[1:4]
    output_tsv = sys.argv[4]
    output_heatmap = sys.argv[5]
    
    main(input_files, output_tsv, output_heatmap)