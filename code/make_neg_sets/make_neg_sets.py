import argparse
import pandas as pd
import random
from Levenshtein import distance as levenshtein_distance


def precompute_allowed_mirnas(positive_file_path, min_allowed_distance):
    # Read only the necessary columns into a pandas DataFrame
    df = pd.read_csv(positive_file_path, delimiter="\t", usecols=['seq.m', 'noncodingRNA_fam'])
    
    # Extract unique pairs of 'seq.m' and 'noncodingRNA_fam'
    unique_mirnas_fams = df.drop_duplicates().values.tolist()
    
    # Create a dictionary to store results
    allowed_mirnas = {}
    
    # Iterate over each unique 'seq.m' to compute distances
    for i, (seq_m, noncodingRNA_fam) in enumerate(unique_mirnas_fams):
        per_mirna_allowed_mirnas = []
        for j, (other_seq_m, other_noncodingRNA_fam) in enumerate(unique_mirnas_fams):
            if i != j:  # Skip comparison with itself
                dist = levenshtein_distance(seq_m, other_seq_m)
                if dist > min_allowed_distance:
                    per_mirna_allowed_mirnas.append((other_seq_m, other_noncodingRNA_fam))
        allowed_mirnas[seq_m] = per_mirna_allowed_mirnas
    
    return allowed_mirnas


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


def generate_negative_samples(block, unique_mirnas_fams, num_negatives, min_required_edit_distance):

    negative_samples = []
    negative_blacklist = []
    unsuccessful = 0

    pos_mirnas = []
    for _, row in block.iterrows():
        pos_mirnas.append(block['seq.m']) 

    for i, pos_row in block.iterrows():
        gene = block['seq.g']
        current_pos_mirna = block['seq.m']
        # The number of negative samples to generate for this positive sample is the sum of the number of negative samples to generate and the number of remaining negative samples to generate due to a failed attempt in the previous positive sample
        n = num_negatives + unsuccessful
        for j in range(1, n + 1): # Iterating over 1 to n+1 to easily account for unsuccessful attempts later
            min_edit_distance = 0
            tries = 0
            while min_edit_distance < min_required_edit_distance and tries < 200:
                random_mirna, random_fam = map(str, random.choice(unique_mirnas_fams[current_pos_mirna]))  
                # Check if the random mirna is not already binding to the gene
                if random_mirna in pos_mirnas:
                    continue
                # Check if the random mirna is not already in the blacklist
                if random_mirna in negative_blacklist:
                    continue
                # Get a minimum edit distance from all the positive mirnas of given gene
                min_edit_distance = min([levenshtein_distance(random_mirna, current_pos_mirna) for pos_mirna in pos_mirnas])
                tries += 1
            if tries == 200:
                # Update `unsuccessful` to add the remaining number of negative samples to generate for this positive sample to the next positive sample to maintain the positive to negative ratio
                unsuccessful = n - j
                if i == len(block.index) - 1:
                    print(f"Warning: Failed to generate all negative samples. Missing {unsuccessful} negative samples.")
                break
                
            negative_sample = pos_row.copy()
            negative_sample['label'] = 0 
            negative_sample['seq.m'] = random_mirna
            negative_sample['noncodingRNA_fam'] = random_fam
            negative_samples.append(negative_sample)

            negative_blacklist.append(random_mirna)
            unsuccessful = 0
    return pd.DataFrame(negative_samples)


def main():
    parser = argparse.ArgumentParser(description="Generate negative samples with specific edit distance.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file name")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    parser.add_argument('--neg_ratio', type=int, default=100, help="Number of negative samples to generate per positive sample")
    parser.add_argument('--min_required_edit_distance', type=int, default=3, help="Minimum required edit distance for negative samples")
    args = parser.parse_args()

    # Set a fixed random seed for reproducibility
    random.seed(42)

    # Load positive samples from the input file
    positive_samples = pd.read_csv(args.ifile, sep='\t')

    # Group by 'gene' and count the number of rows per gene
    gene_count = positive_samples.groupby('seq.g').size().reset_index(name='count')

    # Merge the count information back with the original DataFrame
    positive_samples_gene_count = positive_samples.merge(gene_count, on='seq.g')

    # Sort the DataFrame by the 'count' column in ascending order
    sorted_positive_samples = positive_samples_gene_count.sort_values(by='count', ascending=True)

    with open(sorted_positive_samples) as file_handler:
        unique_mirnas_fams = precompute_allowed_mirnas(sorted_positive_samples, args.min_required_edit_distance)
        for block in yield_gene_blocks(file_handler):
            negatives = generate_negative_samples(block, unique_mirnas_fams, args.neg_ratio, args.min_required_edit_distance)
            negatives.to_csv(args.ofile, sep='\t', mode='a', index=False, header=False)
    
    positive_samples.to_csv(args.ofile, sep='\t', mode='a', index=False, header=False)

if __name__ == "__main__":
    main()

