import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def calculate_percentages(df):
    return (df.div(df.sum()) * 100).round(2)

def plot_mirna_families(df, output_file):
    # get first 10 families from the input
    top_10_families = df.index[:10].tolist()

    # set up the color palette
    colors = ['#f0e442', '#0072b2', '#cc79a7']

    # create the plot with adjusted figure size
    fig, ax = plt.subplots(figsize=(30, 10), dpi=300)

    x = np.arange(len(top_10_families))
    width = 0.25

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate miRNA family plot from TSV data.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file containing miRNA family counts")
    parser.add_argument("-o", "--output", required=True, help="Output PNG file for the generated plot")
    parser.add_argument("-t", "--tsv", required=True, help="Output TSV file for the plot results")
    return parser.parse_args()

def read_input_data(input_file):
    return pd.read_csv(input_file, sep='\t', index_col='miRNA Family')

def prepare_output_data(df, df_percentages, top_10_families):
    output_df = pd.DataFrame(index=top_10_families)
    for column in df.columns:
        output_df[f'{column}_Count'] = df.loc[top_10_families, column]
        output_df[f'{column}_Percentage'] = df_percentages.loc[top_10_families, column]
    return output_df.round(2)

def save_output_data(output_df, output_file):
    output_df.to_csv(output_file, sep='\t')
    print(f"Plot results saved as {output_file}")

  # plot bars centered around the tick
    def plot_centered_bars(ax, x, df, width, colors):
        x_centers = []
        for i, family in enumerate(top_10_families):
            values = df.loc[family]
            non_zero = values[values > 0]
            n_bars = len(non_zero)

            if n_bars == 0:
                x_centers.append(x[i])
                continue

            total_width = n_bars * width
            start = x[i] - total_width / 2 + width / 2

            center = start + (total_width - width) / 2
            x_centers.append(center)

            for j, (col, val) in enumerate(non_zero.items()):
                ax.bar(start, val, width, label=col if i == 0 else "", 
                       color=colors[j % len(colors)], edgecolor='black', linewidth=0.5)
                start += width

        return x_centers

    x_centers = plot_centered_bars(ax, x, df.loc[top_10_families], width, colors)

    # remove labels and title
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_title('')

    # set x-ticks to the center of each group
    ax.set_xticks(x_centers)
    ax.set_xticklabels(top_10_families, rotation=0, ha='center', fontsize=34, weight='semibold')

    # adjust y-axis to 50% max
    ax.set_ylim(0, 50)
    ax.set_yticks(range(0, 51, 10))
    ax.set_yticklabels([f'{i}%' for i in range(0, 51, 10)], fontsize=34)
    ax.yaxis.grid(True, linestyle=':', alpha=0.8, color='black')

    # increase width of x and y axis lines
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    # remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # move x-axis ticks down
    ax.xaxis.set_tick_params(pad=10)

    # add legend without border and with semibold text
    legend = ax.legend(fontsize=34, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)

    # make legend text semibold
    for text in legend.get_texts():
        text.set_fontweight('semibold')

    plt.tight_layout()

    # save the plot as a PNG file
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")

    return top_10_families

def main():
    args = parse_arguments()
    df = read_input_data(args.input)
    df_percentages = calculate_percentages(df)
    top_10_families = plot_mirna_families(df_percentages, args.output)
    output_df = prepare_output_data(df, df_percentages, top_10_families)
    save_output_data(output_df, args.tsv)

if __name__ == "__main__":
    main()