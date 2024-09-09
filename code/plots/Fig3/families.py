import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from collections import OrderedDict
from seed_utils import find_seed_match, reverse_complement

def process_datasets(datasets):
    all_data = {}
    for name, file_path in datasets.items():
        df = pd.read_csv(file_path, sep='\t')
        df['seed_type'] = df.apply(lambda row: find_seed_match(row['seq.g'], row['seq.m']), axis=1)
        all_data[name] = df

    combined_df = pd.concat(all_data.values())

    # Calculate total counts per family across all datasets
    family_counts = combined_df['miRNA_fam'].value_counts()
    total_percentage = (family_counts / family_counts.sum()) * 100

    # Create a DataFrame for total counts and total percentages
    total_counts_df = pd.DataFrame({
        'miRNA_fam': family_counts.index,
        'Total_Count': family_counts.values,
        'Total_Percentage': total_percentage.values
    })

    # Calculate counts and percentages for each dataset
    for dataset_name, df in all_data.items():
        counts = df['miRNA_fam'].value_counts()
        percentages = df['miRNA_fam'].value_counts(normalize=True) * 100

        total_counts_df[f'{dataset_name}_Count'] = total_counts_df['miRNA_fam'].map(counts).fillna(0).astype(int)
        total_counts_df[f'{dataset_name}_Percentage'] = total_counts_df['miRNA_fam'].map(percentages).fillna(0)

    # Sort the DataFrame by total counts in descending order
    total_counts_df = total_counts_df.sort_values(by='Total_Count', ascending=False)

    return combined_df, total_counts_df

def calculate_seed_percentages(df, top_families):
    df_top = df[df['miRNA_fam'].isin(top_families)]

    # Calculate seed type percentages for the top families
    seed_percentages = df_top.groupby('miRNA_fam')['seed_type'].value_counts(normalize=True).unstack() * 100

    # Calculate total seed percentage (excluding 'none')
    seed_percentages['TotalCanonicalSeed'] = seed_percentages['Seed6mer'] + seed_percentages['Seed7mer'] + seed_percentages['Seed8mer']

    # Reorder columns
    seed_percentages = seed_percentages[['Seed6mer', 'Seed7mer', 'Seed8mer', 'SeedNonCanonical', 'TotalCanonicalSeed']]

    # Sort the DataFrame based on the order of top_families
    seed_percentages = seed_percentages.reindex(top_families)

    return seed_percentages

def plot_seed_prevalence(seed_percentages, output_file):
    # Set the desired dimensions in pixels and DPI
    desired_width_px = 7000
    desired_height_px = 3000
    dpi = 300

    # Calculate the figure size in inches
    width_inches = desired_width_px / dpi
    height_inches = desired_height_px / dpi

    # Calculate scaling factor based on the width
    scale_factor = width_inches / 104

    fig, ax = plt.subplots(figsize=(width_inches, height_inches))

    colors = ['#f0e442', '#009e73', '#cc79a7', '#0072b2']

    # Increase bar width even more
    bar_width = 0.0693 * scale_factor * 3.5  # Increased from 2.5 to 3.5

    # Add small space between bars within a group to prevent overlap
    bar_spacing = bar_width * 0.03  # Reduced spacing percentage to allow for wider bars

    group_width = 4 * (bar_width + bar_spacing) - bar_spacing
    # Adjust spacing between grouped columns
    group_spacing = 0.3 * scale_factor * (7000 / 8000)  # Slightly reduced group spacing

    # Reduce padding between y-axis and first group
    left_padding = 0.01 * scale_factor

    index = [left_padding + i * (group_width + group_spacing) for i in range(len(seed_percentages.index))]

    for i, seed_type in enumerate(['SeedNonCanonical', 'Seed6mer', 'Seed7mer', 'Seed8mer']):
        values = seed_percentages[seed_type]
        bar_positions = [x + i*(bar_width + bar_spacing) for x in index]
        bars = ax.bar(bar_positions, values, bar_width, color=colors[i], edgecolor='black', linewidth=0.5*scale_factor, label=seed_type)

        for bar in bars:
            x = bar.get_x()
            y = bar.get_y()
            height = bar.get_height()
            # Increase border width for individual columns
            ax.add_patch(Rectangle((x, y), bar_width, height, fill=False, edgecolor='black', linewidth=4*scale_factor))

    # Bolder horizontal lines for total seed percentage and vertical dotted lines
    for i, total in enumerate(seed_percentages['TotalCanonicalSeed']):
        group_center = index[i] + group_width / 2
        ax.plot([index[i], index[i] + group_width], [total, total], color='black', linewidth=12*scale_factor)
        # Add vertical dotted line
        ax.plot([group_center, group_center], [0, total], color='black', linestyle=':', linewidth=4*scale_factor)

    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Adjust x-ticks to match the center of the groups
    ax.set_xticks([x + (group_width / 2) for x in index])
    ax.set_xticklabels(seed_percentages.index, rotation=0, ha='center', fontsize=112*scale_factor, weight='bold')

    ax.set_ylim(0, 60)
    ax.set_yticks(range(0, 61, 10))
    ax.set_yticklabels([f'{i}%' for i in range(0, 61, 10)], fontsize=112*scale_factor)
    ax.yaxis.grid(True, linestyle=':', alpha=0.8, color='black')

    # Make x and y axes black and bolder
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_linewidth(6*scale_factor)
    ax.spines['bottom'].set_linewidth(6*scale_factor)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x', colors='black', width=6*scale_factor, length=10*scale_factor)
    ax.tick_params(axis='y', colors='black', width=6*scale_factor, length=10*scale_factor)

    ax.xaxis.set_tick_params(pad=15*scale_factor)

    square_size = 96 * scale_factor
    legend_handles = [Patch(facecolor=color, edgecolor='black') for color in colors]

    # Add black line handle for TotalCanonicalSeed
    total_seed_handle = Line2D([0], [0], color='black', linewidth=12*scale_factor)
    legend_handles.append(total_seed_handle)

    legend_labels = ['SeedNonCanonical', 'Seed6mer', 'Seed7mer', 'Seed8mer', 'TotalCanonicalSeed']
    legend = ax.legend(legend_handles, legend_labels, 
                       loc='upper center', bbox_to_anchor=(0.5, -0.05),
                       ncol=5, columnspacing=0.5, handlelength=1, handleheight=1,
                       prop={'size': 96*scale_factor, 'weight': 'semibold'})

    for handle in legend.get_patches():
        handle.set_height(square_size)
        handle.set_width(square_size)

    for text, handle in zip(legend.get_texts(), legend.get_patches()):
        text.set_ha('center')
        text.set_position((handle.get_x() + handle.get_width()/2, text.get_position()[1]))

    legend.get_frame().set_edgecolor('none')

    plt.tight_layout()

    # Ensure the figure size is correct before saving
    fig = plt.gcf()
    fig.set_size_inches(width_inches, height_inches)

    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    print(f"Plot saved as {output_file}")
    print(f"Image dimensions: {desired_width_px}x{desired_height_px} pixels")

def main():
    parser = argparse.ArgumentParser(description='Generate miRNA seed prevalence plot for three datasets.')
    parser.add_argument('dataset1', help='Path to first dataset')
    parser.add_argument('dataset2', help='Path to second dataset')
    parser.add_argument('dataset3', help='Path to third dataset')
    parser.add_argument('-o', '--output', default='seed_prevalence_plot.png', help='Output file name for plot (default: seed_prevalence_plot.png)')
    parser.add_argument('-t', '--tsv_output', default='miRNA_family_analysis.tsv', help='Output TSV file name (default: miRNA_family_analysis.tsv)')

    args = parser.parse_args()

    datasets = OrderedDict([
        (args.dataset1, args.dataset1),
        (args.dataset2, args.dataset2),
        (args.dataset3, args.dataset3)
    ])

    # Process datasets to generate the combined DataFrame with counts and percentages
    combined_df, total_counts_df = process_datasets(datasets)

    # Save the combined results to a TSV file
    total_counts_df.to_csv(args.tsv_output, sep='\t', index=False)
    print(f"Combined TSV file saved as {args.tsv_output}")

    # Generate a plot for the top 10 families based on total counts
    top_10_families = total_counts_df['miRNA_fam'].head(10).tolist()
    seed_percentages = calculate_seed_percentages(combined_df, top_10_families)
    plot_seed_prevalence(seed_percentages, args.output)

if __name__ == "__main__":
    main()
