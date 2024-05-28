import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def plot_data(input_file, plot_file, tsv_file):
    # load the data
    if input_file:
        data = pd.read_csv(input_file, sep='\t')
    else:
        data = pd.read_csv(sys.stdin, sep='\t')

    # extracting relevant columns for the plot
    x = data['mean_total_interactions']
    y = data['mean_percentage']

    # transform x to logarithmic base 10 for better spacing
    x_log10 = np.log10(x + 1)  # adding 1 to avoid log(0)

    # creating the plot
    plt.figure(figsize=(10, 6))
    plt.scatter(x_log10, y)

    # adding labels and title
    plt.xlabel('Log10(Number of Occurrences of each miRNA)')
    plt.ylabel('Percentage of Perfect Seed Interaction')
    plt.title('Number of Interactions vs. Percentage of Perfect Seed Interaction')

    # customize x-axis ticks to show original scale for clarity
    ticks = np.array([1, 10, 100, 1000, 10000])
    tick_labels = [str(tick) for tick in ticks]
    plt.xticks(np.log10(ticks + 1), tick_labels)

    # adding grid lines maybe for better vis
    plt.grid(True, which="both", ls="--")

    # save the plot
    if plot_file:
        plt.savefig(plot_file)
    plt.close()
    
    # save the plot data to a TSV file sorted by 'mean_total_interactions' in descending order
    plot_data = pd.DataFrame({
        'log_mean_total_interactions': x_log10,
        'mean_percentage': y,
        'mean_total_interactions': x
    }).sort_values(by='mean_total_interactions', ascending=False)
    
    if tsv_file:
        plot_data.to_csv(tsv_file, sep='\t', index=False)
    else:
        plot_data.to_csv(sys.stdout, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Scatter Plot Creator")
    parser.add_argument('--ifile', help="Input TSV file (default: stdin)")
    parser.add_argument('--pfile', help="Output plot image file (default: none)", default=None)
    parser.add_argument('--ofile', help="Output TSV file (default: stdout)")
    args = parser.parse_args()

    plot_data(args.ifile, args.pfile, args.ofile)

if __name__ == "__main__":
    main()
