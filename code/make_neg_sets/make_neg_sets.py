import argparse
import pandas as pd
import random
from Levenshtein import distance as levenshtein_distance
import sys
import time

# Get allowed mirnas i.e. mirnas that have an edit distance greater than min_allowed_distance from the given mirna
def precompute_allowed_mirnas(positive_samples_df, min_allowed_distance):
    unique_mirnas = positive_samples_df['noncodingRNA'].unique()
    allowed_mirnas = {}
    
    for mirna in unique_mirnas:
        allowed_mirnas[mirna] = [other_mirna for other_mirna in unique_mirnas 
                                 if mirna != other_mirna and 
                                 levenshtein_distance(mirna, other_mirna) > min_allowed_distance]
    
    return allowed_mirnas

# Get unique pairs of noncodingRNA sequence (seqm) and family
def get_unique_seqm_fam_pairs(positive_samples_df):
    unique_seqm_fam_pairs = positive_samples_df[['noncodingRNA', 'noncodingRNA_fam']].drop_duplicates()
    unique_seqm_fam_pairs_dict = unique_seqm_fam_pairs.set_index('noncodingRNA')['noncodingRNA_fam'].to_dict()
    return unique_seqm_fam_pairs_dict

# Yield blocks of positive examples with the same gene to process at a time, for memory efficiency
def yield_gene_blocks(positive_file_path):
    current_block = []
    current_seq_g = None

    with open(positive_file_path, 'r') as file:
        header_columns = file.readline().strip().split('\t')
        seq_g_index = header_columns.index('gene')

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

# Generate negative examples for a gene block of positive examples
def generate_negative_samples(block, neg_ratio, unique_seqm_fam_pairs_dict, allowed_mirnas, unsuccessful):
    
    # Check if 'feature' and 'test' are consistent within the block
    if block['feature'].nunique() != 1 or block['test'].nunique() != 1:
        return [], unsuccessful
    
    neg_label = 0
    gene = block['gene'].iloc[0]
     
     # Get the set of miRNAs that share this gene as positive examples
    pos_mirnas = block['noncodingRNA'].unique().tolist()
    
    # Get the set of miRNAs that are allowed to be paired with this gene as negative examples
    gene_allowed_mirnas = set.intersection(*[set(allowed_mirnas[mirna]) for mirna in pos_mirnas])
    gene_allowed_mirnas = list(gene_allowed_mirnas)

    # If neg_ratio is 'max', use all the negative examples possible for this gene
    if neg_ratio == 'max':
        unsuccessful = 0
        n_negative_mirnas = gene_allowed_mirnas
    # If neg_ratio is an integer, compute the number of negative examples to generate for this gene based on the neg_ratio, the number of positive examples that share this gene, and the number of unsuccessful negative examples from the previous block. 
    else:
        neg_ratio = int(neg_ratio)
        n = neg_ratio * block.shape[0] + unsuccessful
        # If there are not enough allowed mirnas for this gene, use all of them and record the number of unsuccessful negative examples
        if n > len(gene_allowed_mirnas):
            unsuccessful = n - len(gene_allowed_mirnas)
            n_negative_mirnas = gene_allowed_mirnas
        # If there are enough allowed mirnas for this gene, randomly example n negative examples from the allowed mirnas for this gene
        else:
            unsuccessful = 0
            n_negative_mirnas = random.sample(gene_allowed_mirnas, n)

    feature = block['feature'].iloc[0]
    test = block['test'].iloc[0]

    # Construct the df rows for the negative examples
    negative_sample_rows = [
        [gene, neg_mirna, unique_seqm_fam_pairs_dict[neg_mirna], feature, test, neg_label]
        for neg_mirna in n_negative_mirnas
    ]

    return negative_sample_rows, unsuccessful

def main():
    # Record start time
    start = time.time()

    parser = argparse.ArgumentParser(description="Generate negative examples with specific edit distance.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file name, must be sorted by 'gene'")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    parser.add_argument('--neg_ratio', default=1, help="Number of negative examples to generate per positive example")
    parser.add_argument('--min_required_edit_distance', type=int, default=3, help="Minimum required edit distance for negative examples")
    args = parser.parse_args()

    # Set a fixed random seed for reproducibility
    random.seed(42)  
    
    # Read the entire positive examples file
    positive_samples = pd.read_csv(args.ifile, sep='\t')
    
    # Write header to the output file
    positive_samples.head(0).to_csv(args.ofile, sep='\t', index=False, mode='w')

    unique_seqm_fam_pairs_dict = get_unique_seqm_fam_pairs(positive_samples)
    allowed_mirnas = precompute_allowed_mirnas(positive_samples, args.min_required_edit_distance)

    # Delete the positive_samples dataframe to free up memory
    del positive_samples

    unsuccessful = 0 

    inconsistent_blocks = pd.DataFrame(columns=['gene', 'noncodingRNA', 'noncodingRNA_fam', 'feature', 'test', 'label'])

    with open(args.ofile, 'a') as ofile:
        
        for block in yield_gene_blocks(args.ifile):

            negative_sample_rows, unsuccessful = generate_negative_samples(block, args.neg_ratio, unique_seqm_fam_pairs_dict, allowed_mirnas, unsuccessful)

            if negative_sample_rows:
                # Append positive examples for this block to the output file
                block.to_csv(ofile, sep='\t', index=False, header=False, mode='a')
                # Append negative examples for this block to the output file
                for sublist in negative_sample_rows:
                    ofile.write('\t'.join(map(str, sublist)) + '\n')
            else:
                inconsistent_blocks = pd.concat([inconsistent_blocks, block], ignore_index=True)

    if inconsistent_blocks.shape[0] > 0:
        # Print excluded positive examples block to stderr if no negative examples were generated
        sys.stderr.write(f"Warning: Could not generate negative examples for the following positive examples due to inconsistent feature or chr.g for the same gene. Excluding positive examples. \n")
        sys.stderr.write(inconsistent_blocks.to_string(index=False) + '\n')

    if unsuccessful > 0:
        print(f"Warning: Could not generate {args.neg_ratio} negative examples, missing {unsuccessful} negative examples.")
    
    # Record end time
    end = time.time() 

    # Print the difference between start and end time in secs
    print(f"The time of execution for neg_ratio {args.neg_ratio} is:",(end-start),"s")


if __name__ == "__main__":
    main()
