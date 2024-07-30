import argparse
import pandas as pd
import random
from Levenshtein import distance as levenshtein_distance
import sys

def precompute_allowed_mirnas(positive_samples_df, min_allowed_distance):
    unique_mirnas = positive_samples_df['seq.m'].unique()
    allowed_mirnas = {}
    
    for mirna in unique_mirnas:
        allowed_mirnas[mirna] = [other_mirna for other_mirna in unique_mirnas 
                                 if mirna != other_mirna and 
                                 levenshtein_distance(mirna, other_mirna) > min_allowed_distance]
    
    return allowed_mirnas

def get_unique_seqm_fam_pairs(positive_samples_df):
    unique_seqm_fam_pairs = positive_samples_df[['seq.m', 'noncodingRNA_fam']].drop_duplicates()
    unique_seqm_fam_pairs_dict = unique_seqm_fam_pairs.set_index('seq.m')['noncodingRNA_fam'].to_dict()
    return unique_seqm_fam_pairs_dict

def yield_gene_blocks(positive_file_path):
    current_block = []
    current_seq_g = None

    with open(positive_file_path, 'r') as file:
        header_columns = file.readline().strip().split()
        seq_g_index = header_columns.index('seq.g')

        for line in file:
            columns = line.strip().split()
            seq_g = columns[seq_g_index]

            if seq_g != current_seq_g:
                if current_block:
                    block_df = pd.DataFrame(current_block, columns=header_columns)
                    yield block_df
                current_block = [columns]
                current_seq_g = seq_g
            else:
                current_block.append(columns)

    if current_block:
        block_df = pd.DataFrame(current_block, columns=header_columns)
        yield block_df

def generate_negative_samples(block, num_negatives, unique_seqm_fam_pairs_dict, allowed_mirnas, unsuccessful):
    neg_label = 0
    gene = block['seq.g'].iloc[0]

    pos_mirnas = block['seq.m'].unique().tolist()
    
    # Select allowed miRNAs for positive miRNAs of this gene and intersect them
    gene_allowed_mirnas = set.intersection(*[set(allowed_mirnas[mirna]) for mirna in pos_mirnas])
    gene_allowed_mirnas = gene_allowed_mirnas - set(pos_mirnas)  # Remove positive miRNAs
    gene_allowed_mirnas = list(gene_allowed_mirnas)

    n = num_negatives * block.shape[0] + unsuccessful

    if n > len(gene_allowed_mirnas):
        unsuccessful = n - len(gene_allowed_mirnas)
        n_negative_mirnas = gene_allowed_mirnas
    else:
        unsuccessful = 0
        n_negative_mirnas = random.sample(gene_allowed_mirnas, n)

    # Check if 'feature' and 'test' are consistent within the block
    if block['feature'].nunique() != 1 or block['test'].nunique() != 1:
        print(f"Warning: Inconsistent 'feature' or 'test' values in block for gene {gene}.")
        return [], unsuccessful

    feature = block['feature'].iloc[0]
    test = block['test'].iloc[0]

    negative_sample_rows = [
        [gene, neg_mirna, unique_seqm_fam_pairs_dict[neg_mirna], feature, test, neg_label]
        for neg_mirna in n_negative_mirnas
    ]

    return negative_sample_rows, unsuccessful

def main():
    parser = argparse.ArgumentParser(description="Generate negative samples with specific edit distance.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file name, must be sorted by 'seq.g'")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    parser.add_argument('--neg_ratio', type=int, default=100, help="Number of negative samples to generate per positive sample")
    parser.add_argument('--min_required_edit_distance', type=int, default=3, help="Minimum required edit distance for negative samples")
    args = parser.parse_args()

    random.seed(42)  # Set a fixed random seed for reproducibility
    
    positive_samples = pd.read_csv(args.ifile, sep='\t')
    unique_seqm_fam_pairs_dict = get_unique_seqm_fam_pairs(positive_samples)
    allowed_mirnas = precompute_allowed_mirnas(positive_samples, args.min_required_edit_distance)
    
    negatives_rows = []
    unsuccessful = 0
    for block in yield_gene_blocks(args.ifile):
        negatives, unsuccessful = generate_negative_samples(block, args.neg_ratio, unique_seqm_fam_pairs_dict, allowed_mirnas, unsuccessful)
        negatives_rows.extend(negatives)

    negatives_df = pd.DataFrame(negatives_rows, columns=positive_samples.columns)
    combined_df = pd.concat([positive_samples, negatives_df], ignore_index=True)
    combined_df.to_csv(args.ofile, sep='\t', index=False)

if __name__ == "__main__":
    main()
