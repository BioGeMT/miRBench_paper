import pandas as pd
import sys
import argparse
import matplotlib.pyplot as plt
from collections import Counter

def calculate_feature_prevalence(data):

    # fill nan values in the 'feature' column with an empty string
    data['feature'] = data['feature'].fillna('')

    # split the features into a list of lists
    features_list = data['feature'].apply(lambda x: x.split(',') if x else [])

    # flatten the list of lists into a single list
    all_features = [feature for sublist in features_list for feature in sublist]

    # count the occurrences of each feature
    feature_counts = Counter(all_features)

    # calculate the total number of entries
    total_entries = len(data)

    # create a dataframe to store the results
    feature_prevalence = pd.DataFrame(feature_counts.items(), columns=['Feature', 'Count'])

    # calculate the percentage prevalence
    feature_prevalence['Percentage'] = (feature_prevalence['Count'] / total_entries) * 100

    return feature_prevalence

def read_input(input_file):
    if input_file:
        data = pd.read_csv(input_file, sep='\t')
    else:
        data = pd.read_csv(sys.stdin, sep='\t')
    return data

def write_output(data, output_file):
    if output_file:
        data.to_csv(output_file, sep='\t', index=False)
    else:
        data.to_csv(sys.stdout, sep='\t', index=False)

def plot_feature_prevalence(data, plot_file):
    plt.figure(figsize=(10, 6))
    bars = plt.bar(data['Feature'], data['Count'], color='blue', alpha=0.5)
    plt.xlabel('Feature')
    plt.ylabel('Count')
    plt.title('Feature Prevalence')
    
    for bar, percentage in zip(bars, data['Percentage']):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() / 2, f'{percentage:.2f}%', ha='center', va='center', color='black')
    
    plt.tight_layout()
    
    if plot_file:
        plt.savefig(plot_file)

def main():
    parser = argparse.ArgumentParser(description="feature prevalence calculator")
    parser.add_argument('--ifile', help="input file (default: stdin)")
    parser.add_argument('--ofile', help="output file (default: stdout)")
    parser.add_argument('--pfile', help="plot file (default: none)", default=None)

    args = parser.parse_args()

    # read input data
    data = read_input(args.ifile)

    # process data
    feature_prevalence = calculate_feature_prevalence(data)

    # write output data
    write_output(feature_prevalence, args.ofile)

    # plot feature prevalence
    plot_feature_prevalence(feature_prevalence, args.pfile)

if __name__ == "__main__":
    main()
