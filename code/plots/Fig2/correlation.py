import sys
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def load_data(file_path):
    return pd.read_csv(file_path, sep='\t')

def calculate_mirna_percentages(df):
    mirna_counts = df['noncodingRNA_fam'].value_counts()
    total_count = mirna_counts.sum()
    return (mirna_counts / total_count) * 100

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
    # set the figure size
    fig_size_inches = 2800 / 300

    # create a new figure with the specified size
    fig = plt.figure(figsize=(fig_size_inches, fig_size_inches))
    
    # create a grid specification for the figure, with 1 row and 22 columns
    gs = fig.add_gridspec(1, 22)
    
    # add a subplot for the heatmap, using the first 19 columns
    ax_heatmap = fig.add_subplot(gs[0, :19])

    # define the color scheme for the heatmap
    colors = ['#f0e442', '#d55e00']
    cmap = LinearSegmentedColormap.from_list("custom", colors)
    
    # create a mask to show only the upper triangle of the correlation matrix
    mask = np.tril(np.ones_like(corr_matrix, dtype=bool), k=-1)
    
    # get the minimum and maximum correlation values for color scaling
    vmin = 0.5 #corr_matrix.min().min()
    vmax = corr_matrix.max().max()

    # create the heatmap 
    heatmap = sns.heatmap(corr_matrix, ax=ax_heatmap, mask=mask, annot=True, fmt=".3f", cmap=cmap, 
                          vmin=vmin, vmax=vmax, square=True, cbar=False,
                          annot_kws={'size': 22, 'weight': 'medium'})

    # add black borders to each cell in the heatmap
    for i in range(len(corr_matrix)):
        for j in range(len(corr_matrix)):
            if not mask[i, j]:
                heatmap.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='none', lw=2))

    # set the rotation for x-axis and y-axis labels
    ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), rotation=0)
    ax_heatmap.set_yticklabels(ax_heatmap.get_yticklabels(), rotation=90)

    # set the font size and weight for tick labels
    tick_font_size = 18
    ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), weight='medium', size=tick_font_size)
    ax_heatmap.set_yticklabels(ax_heatmap.get_yticklabels(), weight='medium', size=tick_font_size)

    # adjust the padding for tick labels
    ax_heatmap.tick_params(axis='x', which='major', pad=7)
    ax_heatmap.tick_params(axis='y', which='major', pad=7)

    # set the linewidth for the heatmap border
    for spine in ax_heatmap.spines.values():
        spine.set_linewidth(1.5)

    # add a subplot for the colorbar, using the 21st column
    ax_colorbar = fig.add_subplot(gs[0, 20])

    # calculate the dimensions for the colorbar
    heatmap_height = ax_heatmap.get_position().height
    box_height = heatmap_height / len(corr_matrix)
    colorbar_height = 3 * box_height
    colorbar_bottom = 0.5 - (colorbar_height / 2)

    # set the position and size of the colorbar
    colorbar_width = ax_colorbar.get_position().width * 0.75
    ax_colorbar.set_position([ax_colorbar.get_position().x0,
                              colorbar_bottom,
                              colorbar_width,
                              colorbar_height])

    # create a ScalarMappable for the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    
    # add the colorbar to the figure
    cbar = fig.colorbar(sm, cax=ax_colorbar)
    cbar.ax.tick_params(labelsize=20)

    # remove the outline of the colorbar
    cbar.outline.set_visible(False)

    # save the figure as an image file
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    plt.close()

def parse_arguments():
    if len(sys.argv) != 6:
        print("Usage: python script.py input1.tsv input2.tsv input3.tsv output_correlations.tsv output_heatmap.png")
        sys.exit(1)
    return sys.argv[1:4], sys.argv[4], sys.argv[5]

def process_input_files(input_files):
    datasets = [calculate_mirna_percentages(load_data(file)) for file in input_files]
    file_names = [os.path.basename(file) for file in input_files]
    return datasets, file_names

def save_correlation_matrix(corr_matrix, output_tsv):
    corr_matrix.to_csv(output_tsv, sep='\t')

def main():
    input_files, output_tsv, output_heatmap = parse_arguments()
    datasets, file_names = process_input_files(input_files)
    corr_matrix = calculate_correlations(datasets, file_names)
    np.fill_diagonal(corr_matrix.values, 1)
    save_correlation_matrix(corr_matrix, output_tsv)
    create_heatmap(corr_matrix, output_heatmap)

if __name__ == "__main__":
    main()