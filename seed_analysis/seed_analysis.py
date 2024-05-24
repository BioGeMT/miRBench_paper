import os
import pandas as pd
from Bio.Seq import Seq
import sys

def load_tsv_files(input_folder):
    # read all TSV files from the input folder
    dataframes = []
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".tsv"):
            file_path = os.path.join(input_folder, file_name)
            df = pd.read_csv(file_path, sep='\t')
            dataframes.append(df)
    return pd.concat(dataframes, ignore_index=True)

def reverse_complement(seq):
    # return the reverse complement of a given sequence
    return str(Seq(seq).reverse_complement())

def perfect_match(mirna, target):
    # check if the seed region (nucleotides 2 to 7) of the microRNA perfectly matches the reverse complement of the target sequence
    mirna_seed = mirna[1:7]
    target_rc = reverse_complement(target)
    has_match = mirna_seed in target_rc
    return has_match

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

def save_results(aggregate_df, summary_stats, output_folder):
    # save the aggregated results and summary statistics to TSV files
    aggregate_df.to_csv(os.path.join(output_folder, 'aggregate_results.tsv'), sep='\t', index=False)
    summary_stats.to_csv(os.path.join(output_folder, 'summary_stats.tsv'), sep='\t', index=False)

def main(input_folder, output_folder):
    # load all TSV files from the input folder
    dataframes = load_tsv_files(input_folder)

    # process the file to get interaction counts and perfect matches
    results = process_file(dataframes)

    # aggregate the results and calculate summary statistics
    aggregate_df, summary_stats = aggregate_results(results)

    # save the results to the output folder
    save_results(aggregate_df, summary_stats, output_folder)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python seed_analysis.py <input_folder> <output_folder>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    main(input_folder, output_folder)

