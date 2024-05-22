import os
import pandas as pd
from Bio.Seq import Seq
import argparse
import sys

def load_tsv_files(input_stream):
    # Read from the input stream (which could be stdin or files from a folder)
    dataframes = pd.read_csv(input_stream, sep='\t')
    return dataframes

def reverse_complement(seq):
    # return the reverse complement of a given sequence
    return str(Seq(seq).reverse_complement())

def perfect_match(mirna, target):
    # check if the seed region (nucleotides 2 to 7) of the microRNA perfectly matches the reverse complement of the target sequence
    mirna_seed = mirna[1:7]
    target_rc = reverse_complement(target)
    return mirna_seed in target_rc

def process_file(df):
    # apply the perfect_match function to each row and create a new column 'perfect_match'
    df['perfect_match'] = df.apply(lambda row: perfect_match(row['seq.m'], row['seq.g']), axis=1)
    
    # group by noncodingRNA_fam and count total interactions and perfect matches
    result = df.groupby('noncodingRNA_fam').agg(
        total_interactions=('seq.m', 'count'),
        perfect_matches=('perfect_match', 'sum')
    ).reset_index()
    
    # calculate the percentage of perfect matches
    result['percentage'] = (result['perfect_matches'] / result['total_interactions']) * 100
    return result

def aggregate_results(results):
    # aggregate results from multiple dataframes and calculate summary statistics
    aggregate_df = results
    summary_stats = aggregate_df.groupby('noncodingRNA_fam').agg(
        mean_total_interactions=('total_interactions', 'mean'),
        std_total_interactions=('total_interactions', 'std'),
        median_total_interactions=('total_interactions', 'median'),
        min_total_interactions=('total_interactions', 'min'),
        max_total_interactions=('total_interactions', 'max'),
        q1_total_interactions=('total_interactions', lambda x: x.quantile(0.25)),
        q3_total_interactions=('total_interactions', lambda x: x.quantile(0.75)),
        iqr_total_interactions=('total_interactions', lambda x: x.quantile(0.75) - x.quantile(0.25)),
        variance_total_interactions=('total_interactions', 'var'),
        mean_perfect_matches=('perfect_matches', 'mean'),
        std_perfect_matches=('perfect_matches', 'std'),
        median_perfect_matches=('perfect_matches', 'median'),
        min_perfect_matches=('perfect_matches', 'min'),
        max_perfect_matches=('perfect_matches', 'max'),
        q1_perfect_matches=('perfect_matches', lambda x: x.quantile(0.25)),
        q3_perfect_matches=('perfect_matches', lambda x: x.quantile(0.75)),
        iqr_perfect_matches=('perfect_matches', lambda x: x.quantile(0.75) - x.quantile(0.25)),
        variance_perfect_matches=('perfect_matches', 'var'),
        mean_percentage=('percentage', 'mean'),
        std_percentage=('percentage', 'std'),
        median_percentage=('percentage', 'median'),
        min_percentage=('percentage', 'min'),
        max_percentage=('percentage', 'max'),
        q1_percentage=('percentage', lambda x: x.quantile(0.25)),
        q3_percentage=('percentage', lambda x: x.quantile(0.75)),
        iqr_percentage=('percentage', lambda x: x.quantile(0.75) - x.quantile(0.25)),
        variance_percentage=('percentage', 'var')
    ).reset_index()
    return aggregate_df, summary_stats

def save_results(aggregate_df, summary_stats, output_stream):
    # save the aggregated results and summary statistics to TSV files or to stdout
    aggregate_df.to_csv(output_stream, sep='\t', index=False)
    summary_stats.to_csv(output_stream, sep='\t', index=False)

def main():
    # argument parsing
    parser = argparse.ArgumentParser(description="miRNA Seed Analysis")
    parser.add_argument('--ifolder', help="Input folder containing TSV files or STDIN (default: STDIN)", default=None)
    parser.add_argument('--ofolder', help="Output folder for results or STDOUT (default: STDOUT)", default=None)
    args = parser.parse_args()

    # load all TSV files from the input folder or stdin
    if args.ifolder:
        dataframes = load_tsv_files(open(os.path.join(args.ifolder), 'r'))
    else:
        dataframes = load_tsv_files(sys.stdin)

    # process the file to get interaction counts and perfect matches
    results = process_file(dataframes)

    # aggregate the results and calculate summary statistics
    aggregate_df, summary_stats = aggregate_results(results)

    # save the results to the output folder or stdout
    if args.ofolder:
        with open(os.path.join(args.ofolder, 'aggregate_results.tsv'), 'w') as aggregate_output:
            save_results(aggregate_df, summary_stats, aggregate_output)
    else:
        save_results(aggregate_df, summary_stats, sys.stdout)

if __name__ == "__main__":
    main()
