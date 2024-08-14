import sys
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
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
    # Calculate figure size in inches for 2480x2480 pixels at 300 DPI
    fig_size_inches = 2480 / 300  # Equal to 8.27 inches

    fig = plt.figure(figsize=(fig_size_inches, fig_size_inches))
    
    # Create a gridspec for the heatmap and colorbar
    gs = fig.add_gridspec(1, 22)
    
    # Create an axes for the heatmap
    ax_heatmap = fig.add_subplot(gs[0, :19])
    
    # Define the color range
    colors = ['#ffd966', '#d42b1f']
    
    # Create a custom colormap
    cmap = LinearSegmentedColormap.from_list("custom", colors)
    
    # Mask the upper triangle of the heatmap
    mask = np.tril(np.ones_like(corr_matrix, dtype=bool), k=-1)

    # Find the minimum and maximum correlation values
    vmin = corr_matrix.min().min()
    vmax = corr_matrix.max().max()

    # Create the heatmap with a black border
    heatmap = sns.heatmap(corr_matrix, ax=ax_heatmap, mask=mask, annot=True, fmt=".2f", cmap=cmap, 
                          vmin=vmin, vmax=vmax, square=True, cbar=False,
                          annot_kws={'size': 22, 'weight': 'semibold'},
                          linewidths=0, linecolor='black')

    # Manually add black borders only to non-masked cells
    for i in range(len(corr_matrix)):
        for j in range(len(corr_matrix)):
            if not mask[i, j]:
                heatmap.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=1.5))

    # Adjust x and y axis labels
    ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), rotation=0)
    ax_heatmap.set_yticklabels(ax_heatmap.get_yticklabels(), rotation=90)
    
    tick_font_size = 15
    ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), weight='semibold', size=tick_font_size)
    ax_heatmap.set_yticklabels(ax_heatmap.get_yticklabels(), weight='semibold', size=tick_font_size)

    # Adjust tick parameters
    ax_heatmap.tick_params(axis='x', which='major', pad=7)
    ax_heatmap.tick_params(axis='y', which='major', pad=7)

    # Increase the width of axes lines
    for spine in ax_heatmap.spines.values():
        spine.set_linewidth(1.5)

    # Create a separate axes for the colorbar
    ax_colorbar = fig.add_subplot(gs[0, 20])
    
    heatmap_height = ax_heatmap.get_position().height
    box_height = heatmap_height / len(corr_matrix)
    colorbar_height = 3 * box_height
    colorbar_bottom = 0.5 - (colorbar_height / 2)
    
    colorbar_width = ax_colorbar.get_position().width * 0.75
    ax_colorbar.set_position([ax_colorbar.get_position().x0,
                              colorbar_bottom,
                              colorbar_width,
                              colorbar_height])

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=ax_colorbar)
    cbar.ax.tick_params(labelsize=19)
    
    cbar.outline.set_visible(False)

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