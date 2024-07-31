import argparse
import pandas as pd
import random
from Levenshtein import distance as levenshtein_distance
import sys
import time

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
        header_columns = file.readline().strip().split('\t')
        seq_g_index = header_columns.index('seq.g')

        for line in file:
            columns = line.strip().split('\t')
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

def generate_negative_samples(block, neg_ratio, unique_seqm_fam_pairs_dict, allowed_mirnas, unsuccessful):
    neg_label = 0
    gene = block['seq.g'].iloc[0]

    pos_mirnas = block['seq.m'].unique().tolist()
    
    # Get the set of miRNAs that are allowed to be negative samples for this gene
    gene_allowed_mirnas = set.intersection(*[set(allowed_mirnas[mirna]) for mirna in pos_mirnas])
    gene_allowed_mirnas = list(gene_allowed_mirnas)

    n = neg_ratio * block.shape[0] + unsuccessful

    if n > len(gene_allowed_mirnas):
        unsuccessful = n - len(gene_allowed_mirnas)
        n_negative_mirnas = gene_allowed_mirnas
    else:
        unsuccessful = 0
        n_negative_mirnas = random.sample(gene_allowed_mirnas, n)

    # Check if 'feature' and 'test' are consistent within the block
    if block['feature'].nunique() != 1 or block['test'].nunique() != 1:
        print(f"Warning: Inconsistent 'feature' or 'test' values in block for gene {gene}.")
        return [], unsuccessful # or ignore block?

    feature = block['feature'].iloc[0]
    test = block['test'].iloc[0]

    negative_sample_rows = [
        [gene, neg_mirna, unique_seqm_fam_pairs_dict[neg_mirna], feature, test, neg_label]
        for neg_mirna in n_negative_mirnas
    ]

    return negative_sample_rows, unsuccessful

def main():
    # record start time
    start = time.time()

    parser = argparse.ArgumentParser(description="Generate negative samples with specific edit distance.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file name, must be sorted by 'seq.g'")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    parser.add_argument('--neg_ratio', type=int, default=100, help="Number of negative samples to generate per positive sample")
    parser.add_argument('--min_required_edit_distance', type=int, default=3, help="Minimum required edit distance for negative samples")
    args = parser.parse_args()

    random.seed(42)  # Set a fixed random seed for reproducibility
    
    # Read the entire positive samples file
    positive_samples = pd.read_csv(args.ifile, sep='\t')
    
    # Write positive samples to the output file
    positive_samples.to_csv(args.ofile, sep='\t', index=False, mode='w')

    unique_seqm_fam_pairs_dict = get_unique_seqm_fam_pairs(positive_samples)
    allowed_mirnas = precompute_allowed_mirnas(positive_samples, args.min_required_edit_distance)
    del positive_samples

    unsuccessful = 0 

    with open(args.ofile, 'a') as ofile:
        
        for block in yield_gene_blocks(args.ifile):
            negative_sample_rows, unsuccessful = generate_negative_samples(block, args.neg_ratio, unique_seqm_fam_pairs_dict, allowed_mirnas, unsuccessful)
            
            # Append negative samples for this block to the output file
            for sublist in negative_sample_rows:
                ofile.write('\t'.join(map(str, sublist)) + '\n')

    if unsuccessful > 0:
        print(f"Warning: Could not generate {args.neg_ratio} negative samples, missing {unsuccessful} negative samples.")
    
    # record end time
    end = time.time() 

    # print the difference between start and end time in secs
    print(f"The time of execution for neg_ratio {args.neg_ratio} is:",(end-start),"s")


if __name__ == "__main__":
    main()
