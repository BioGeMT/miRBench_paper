import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys

# Use 'Agg' backend for matplotlib
import matplotlib
matplotlib.use('Agg')

def calculate_family_statistics(data, min_interactions):
    # filter families with the specified minimum interactions
    data = data[data['mean_total_interactions'] >= min_interactions]

    # create a dataframe to store the results
    family_stats = data[['noncodingRNA_fam', 'mean_percentage', 'min_percentage', 'max_percentage']]

    # sort the data by mean percentage for plotting (low to high)
    family_stats_sorted_for_plot = family_stats.sort_values(by='mean_percentage')

    # sort the data by mean percentage for TSV output (high to low)
    family_stats_sorted_for_tsv = family_stats.sort_values(by='mean_percentage', ascending=False)

    return family_stats_sorted_for_plot, family_stats_sorted_for_tsv

def read_input(input_file):
    if input_file:
        data = pd.read_csv(input_file, sep='\t')
    else:
        data = pd.read_csv(sys.stdin, sep='\t')
    return data

def write_output(data, output_file):
    # only save necessary columns to the output tsv
    output_data = data[['noncodingRNA_fam', 'mean_percentage']]
    if output_file:
        output_data.to_csv(output_file, sep='\t', index=False)
    else:
        output_data.to_csv(sys.stdout, sep='\t', index=False)

def plot_family_statistics(data, plot_file):
    plt.figure(figsize=(20, 10))

    # custom boxplot without error bars
    for i in range(data.shape[0]):
        plt.plot([i-0.2, i+0.2], [data['min_percentage'].iloc[i]]*2, color='black', linewidth=3, label='min value' if i == 0 else "")
        plt.plot([i-0.2, i+0.2], [data['max_percentage'].iloc[i]]*2, color='black', linewidth=3, label='max value' if i == 0 else "")
        plt.plot([i-0.2, i+0.2], [data['mean_percentage'].iloc[i]]*2, color='red', linewidth=3, label='mean value' if i == 0 else "")

    for i in range(data.shape[0]):
        plt.gca().add_patch(plt.Rectangle((i-0.2, data['min_percentage'].iloc[i]), 0.4, data['max_percentage'].iloc[i] - data['min_percentage'].iloc[i], 
                                          fill=True, color='grey', alpha=0.4))

    plt.xticks(range(data.shape[0]), data['noncodingRNA_fam'], rotation=90)
    plt.title('box plot of mean percentage per family')
    plt.xlabel('noncoding rna family')
    plt.ylabel('mean percentage')
    plt.ylim(0, 100)  # ensure y-axis is from 0 to 100%
    plt.legend()
    plt.tight_layout()
    
    if plot_file:
        plt.savefig(plot_file)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="box plot creator")
    parser.add_argument('--ifile', help="input file (default: stdin)")
    parser.add_argument('--ofile', help="output tsv file (default: stdout)")
    parser.add_argument('--pfile', help="plot file (default: save plot)", default=None)
    parser.add_argument('--min_interactions', type=int, help="minimum number of interactions to filter families (default: 10)", default=10)

    args = parser.parse_args()

    # read input data
    data = read_input(args.ifile)

    # process data
    family_stats_sorted_for_plot, family_stats_sorted_for_tsv = calculate_family_statistics(data, args.min_interactions)

    # write output data
    write_output(family_stats_sorted_for_tsv, args.ofile)

    # plot family statistics
    if args.pfile:
        plot_family_statistics(family_stats_sorted_for_plot, args.pfile)

if __name__ == "__main__":
    main()
