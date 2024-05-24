import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys

def load_data(input_file):
    # load the data from the input file
    return pd.read_csv(input_file, sep='\t')

def create_box_plot(data, output_file, min_interactions):
    # filter families with the specified minimum interactions
    filtered_data = data[data['mean_total_interactions'] >= min_interactions]

    # sort the data by mean percentage
    filtered_data = filtered_data.sort_values(by='mean_percentage')

    plt.figure(figsize=(20, 10))

    # custom boxplot without error bars
    for i in range(filtered_data.shape[0]):
        plt.plot([i-0.2, i+0.2], [filtered_data['min_percentage'].iloc[i]]*2, color='black', linewidth=3, label='Min Value' if i == 0 else "")
        plt.plot([i-0.2, i+0.2], [filtered_data['max_percentage'].iloc[i]]*2, color='black', linewidth=3, label='Max Value' if i == 0 else "")
        plt.plot([i-0.2, i+0.2], [filtered_data['mean_percentage'].iloc[i]]*2, color='red', linewidth=3, label='Mean Value' if i == 0 else "")

    for i in range(filtered_data.shape[0]):
        plt.gca().add_patch(plt.Rectangle((i-0.2, filtered_data['min_percentage'].iloc[i]), 0.4, filtered_data['max_percentage'].iloc[i] - filtered_data['min_percentage'].iloc[i], 
                                          fill=True, color='grey', alpha=0.4))

    # external only error bars
    # plt.errorbar(x=range(filtered_data.shape[0]), y=filtered_data['min_percentage'] - filtered_data['std_percentage'], 
    #              yerr=filtered_data['std_percentage'], fmt='none', capsize=5, color='black', linestyle='none')
    # plt.errorbar(x=range(filtered_data.shape[0]), y=filtered_data['max_percentage'] + filtered_data['std_percentage'], 
    #              yerr=filtered_data['std_percentage'], fmt='none', capsize=5, color='black', linestyle='none')

    plt.xticks(range(filtered_data.shape[0]), filtered_data['noncodingRNA_fam'], rotation=90)
    plt.title('Box Plot of Mean Percentage per Family')
    plt.xlabel('Noncoding RNA Family')
    plt.ylabel('Mean Percentage')
    plt.ylim(0, 100)  # ensure y-axis is from 0 to 100%
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def main():
    # argument parsing
    parser = argparse.ArgumentParser(description="Box Plot Creator")
    parser.add_argument('--ifile', help="Input file (default: STDIN)", default=None)
    parser.add_argument('--ofile', help="Output file (default: STDOUT)", default=None)
    parser.add_argument('--min_interactions', type=int, help="Minimum number of interactions to filter families (default: 10)", default=10)
    args = parser.parse_args()

    # load data from input file or stdin
    if args.ifile:
        data = load_data(args.ifile)
    else:
        data = pd.read_csv(sys.stdin, sep='\t')

    # create box plot and save to output file or stdout
    if args.ofile:
        create_box_plot(data, args.ofile, args.min_interactions)
    else:
        plt.figure()
        create_box_plot(data, 'boxplot.png', args.min_interactions)
        plt.show()

if __name__ == "__main__":
    main()
