import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def plot_data(input_file, output_file):
    # Load the data
    if input_file:
        data = pd.read_csv(input_file, sep='\t')
    else:
        data = pd.read_csv(sys.stdin, sep='\t')

    # Extracting relevant columns for the plot
    x = data['mean_total_interactions']
    y = data['mean_percentage']

    # Transform x to logarithmic base 10 for better spacing
    x_log10 = np.log10(x + 1)  # Adding 1 to avoid log(0)

    # Creating the plot
    plt.figure(figsize=(10, 6))
    plt.scatter(x_log10, y)

    # Adding labels and title
    plt.xlabel('Log10(Number of Occurrences of each miRNA)')
    plt.ylabel('Percentage of Perfect Seed Interaction')
    plt.title('Number of Interactions vs. Percentage of Perfect Seed Interaction')

    # Customize x-axis ticks to show original scale for clarity
    ticks = np.array([1, 10, 100, 1000, 10000,])
    tick_labels = [str(tick) for tick in ticks]
    plt.xticks(np.log10(ticks + 1), tick_labels)

    # Adding grid lines for better readability
    plt.grid(True, which="both", ls="--")

    # Save the plot
    if output_file:
        plt.savefig(output_file)
    else:
        plt.savefig(sys.stdout.buffer, format='png')

def main():
    parser = argparse.ArgumentParser(description='Plot Number of Interactions vs. Percentage of Perfect Seed Interaction')
    parser.add_argument('--ifile', type=str, help='Input TSV file', default=None)
    parser.add_argument('--ofile', type=str, help='Output image file', default=None)
    args = parser.parse_args()

    plot_data(args.ifile, args.ofile)

if __name__ == "__main__":
    main()
