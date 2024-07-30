import argparse
import pandas as pd
import random
from Levenshtein import distance as levenshtein_distance
import sys

def get_unique_seqm_fam_pairs(positive_file_path):

    # Read only the necessary columns into a pandas DataFrame
    df = pd.read_csv(positive_file_path, delimiter="\t", usecols=['seq.m', 'noncodingRNA_fam'])

    # Extract unique pairs of 'seq.m' and 'noncodingRNA_fam'
    unique_seqm_fam_pairs = df.drop_duplicates().values.tolist()

    return unique_seqm_fam_pairs

def yield_gene_blocks(positive_file_path):
    current_block = []
    current_seq_g = None

    with open(positive_file_path, 'r') as file:
        # Read the header (first line)
        header_columns = file.readline().strip().split()
        seq_g_index = header_columns.index('seq.g')

        # Process the rest of the file
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

def compute_allowed_mirnas(gene_positives, min_allowed_distance):

    allowed_mirnas = {}
    
    for mirna in gene_positives:
        per_mirna_allowed_mirnas = []
        for other_mirna in gene_positives:
            if mirna != other_mirna:
                dist = levenshtein_distance(mirna, other_mirna)
                if dist > min_allowed_distance:
                    per_mirna_allowed_mirnas.append(other_mirna)
        allowed_mirnas[mirna] = per_mirna_allowed_mirnas            

    return allowed_mirnas

def intersect_allowed_mirnas(allowed_mirnas):
    return set.intersection(*map(set, allowed_mirnas.values())) # returns a set

def generate_negative_samples(block, num_negatives, unique_seqm_fam_pairs, unsuccessful, min_required_edit_distance):
    
    neg_label = 0
    feature = block['feature']
    test = block['test']
    gene = block['seq.g']

    negative_sample_rows = []

    pos_mirnas = block['seq.m'].unique().tolist()

    gene_allowed_mirnas = compute_allowed_mirnas(pos_mirnas, min_required_edit_distance)

    negative_mirnas = intersect_allowed_mirnas(gene_allowed_mirnas)

    n = num_negatives*block.shape[0] + unsuccessful

    if n > len(negative_mirnas):
        unsuccessful = n - len(negative_mirnas)
    else:
        unsuccessful = 0

    n_negative_mirnas = random.sample(negative_mirnas, n)

    if block[feature].nunique() == 1:
        feature = block[feature].iloc[0]
    else:
        print(f"Warning: Multiple values for 'feature' in block {gene}.") 
        sys.exit(1)

    if block[test].nunique() == 1:
        test = block[test].iloc[0]
    else:
        print(f"Warning: Multiple values for 'test' in block {gene}.")
        sys.exit(1)

    for neg_mirna in n_negative_mirnas:
        neg_row = [gene, neg_mirna, unique_seqm_fam_pairs[neg_mirna], feature, test, neg_label]
        negative_sample_rows.append(neg_row)    

    return negative_sample_rows, unsuccessful


def main():
    parser = argparse.ArgumentParser(description="Generate negative samples with specific edit distance.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file namem, must be sorted by 'seq.g'")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    parser.add_argument('--neg_ratio', type=int, default=100, help="Number of negative samples to generate per positive sample")
    parser.add_argument('--min_required_edit_distance', type=int, default=3, help="Minimum required edit distance for negative samples")
    args = parser.parse_args()

    # Set a fixed random seed for reproducibility
    random.seed(42)

    # Load positive samples from the input file
    positive_samples = pd.read_csv(args.ifile, sep='\t')

    with open(positive_samples) as file_handler:
        unique_seqm_fam_pairs = get_unique_seqm_fam_pairs(file_handler)
        negatives_rows = []
        unsuccessful = 0
        for block in yield_gene_blocks(file_handler):
            negatives, unsuccessful = generate_negative_samples(block, args.neg_ratio, unique_seqm_fam_pairs, unsuccessful, args.min_required_edit_distance)
            negatives_rows.append(negatives)

    negatives_df = pd.DataFrame(negatives_rows)
    combined_df = pd.concat([positive_samples, negatives_df], ignore_index=True)
    combined_df.to_csv(args.ofile, sep='\t', index=False)

if __name__ == "__main__":
    main()

