import argparse
import pandas as pd
import time
import hashlib
from miRBench.encoder import get_encoder
from miRBench.predictor import get_predictor

# Yield blocks of positive examples with the same mirnafam to process at a time
def yield_mirnafam_blocks(positive_file_path):
    current_block = []
    current_mirnafam = None

    with open(positive_file_path, 'r') as file:
        header_columns = file.readline().strip().split('\t')
        mirnafam_index = header_columns.index('noncodingRNA_fam')

        for line in file:
            columns = line.strip().split('\t')
            mirnafam = columns[mirnafam_index]

            if mirnafam != current_mirnafam:
                if current_block:
                    block_df = pd.DataFrame(current_block, columns=header_columns)
                    yield block_df
                current_block = [columns]
                current_mirnafam = mirnafam
            else:
                current_block.append(columns)

    if current_block:
        block_df = pd.DataFrame(current_block, columns=header_columns)
        yield block_df

# Get seed types for the given dataframe
def get_seeds(df):
    seed_types = ["Seed6mer", "Seed6merBulgeOrMismatch"]
    for tool in seed_types:       
        encoder = get_encoder(tool)
        predictor = get_predictor(tool)
        encoded_input = encoder(df)
        output = predictor(encoded_input)
        df[tool] = output
    return df

def make_reproducibility_seed(string)

    # Generate a SHA-256 hash and get the hexadecimal string
    miRNA_hash_hex = hashlib.sha256(string.encode()).hexdigest()
    # Convert the hexadecimal hash to a decimal integer
    miRNA_hash_int = int(miRNA_hash_hex, 16)
    # Reduce the size using modulo (e.g., within the range of a 32-bit unsigned integer) and set it as the seed
    seed = miRNA_hash_int % 4294967295

    return seed

def process_valid_negatives(valid_negatives, block_columns):    

    valid_negatives = valid_negatives.drop(columns=['Seed6mer', 'Seed6merBulgeOrMismatch'])
    valid_negatives['label'] = 0
    valid_negatives = valid_negatives[block_columns]

    return valid_negatives

def process_block(block, positive_samples, all_clusters, output_file, interaction_type):

    # Generate seed for reproducibility
    miRNA_name = block['noncodingRNA_name'].iloc[0]
    seed = make_reproducibility_seed(miRNA_name)

    # Get the set of cluster ids that share this miRNA family
    block_clusters = block['gene_cluster_ID'].unique().tolist()

    # Get the set of cluster ids that are allowed to be paired with this miRNA family
    mirfam_allowed_clusters = [cluster for cluster in all_clusters if cluster not in block_clusters]

    # Pool gene rows from allowed clusters
    negative_gene_pool = positive_samples[positive_samples['gene_cluster_ID'].isin(mirfam_allowed_clusters)]

    # Shuffle the negative gene pool and drop duplicates based on ClusterID
    negative_gene_pool = negative_gene_pool.sample(frac=1, random_state=seed).drop_duplicates(subset=['gene_cluster_ID'], keep='first')
    
    # Make list of unique noncodingRNA values in block
    unique_mirnas = block['noncodingRNA'].unique().tolist()

    # Compute occurrences of each unique miRNA in the block
    mirna_counts = block['noncodingRNA'].value_counts()
    
    # Define columns
    gene_columns = ['gene', 'feature', 'test', 'chr', 'start', 'end', 'strand', 'gene_cluster_ID']
    mirna_columns = ['noncodingRNA', 'noncodingRNA_name', 'noncodingRNA_fam']
    seed_columns = ['Seed6mer', 'Seed6merBulgeOrMismatch']

    for mirna in unique_mirnas:
        # Get frequency of the current miRNA
        mirna_frequency = mirna_counts[mirna]

        # Initialize valid_negatives dataframe
        valid_negatives = pd.DataFrame(columns=gene_columns + mirna_columns + seed_columns)

        # Increment seed for each miRNA
        seed += 1

        # Shuffle the negative gene pool with incrementing seed for each miRNA
        negative_gene_pool = negative_gene_pool.sample(frac=1, random_state=seed)

        # Iterate over each row of the negative gene pool
        for index, row in negative_gene_pool.iterrows():

            # Get gene columns from the row
            negative_candidate = pd.DataFrame([row])[gene_columns].copy()

            # Add miRNA columns to the negative candidate
            negative_candidate['noncodingRNA'] = mirna
            negative_candidate['noncodingRNA_name'] = block[block['noncodingRNA'] == mirna]['noncodingRNA_name'].iloc[0] # Assumes that the name is the same for all occurrences of the miRNA
            negative_candidate['noncodingRNA_fam'] = block[block['noncodingRNA'] == mirna]['noncodingRNA_fam'].iloc[0] # Assumes that the family is the same for all occurrences of the miRNA
            
            # Compute seeds for the negative candidate
            negative_candidate = get_seeds(negative_candidate)

            # Filter the negative candidates based on interaction type
            if interaction_type == 'nonseed':
                negative_candidate = negative_candidate[negative_candidate['Seed6merBulgeOrMismatch'] == 0]
            elif interaction_type == 'canonicalseed':
                negative_candidate = negative_candidate[negative_candidate['Seed6mer'] == 1]
            elif interaction_type == 'noncanonicalseed':
                negative_candidate = negative_candidate[(negative_candidate['Seed6mer'] == 0) & (negative_candidate['Seed6merBulgeOrMismatch'] == 1)]
            
            # If negative candidate is empty
            if negative_candidate.empty:
                continue
            # If negative candidate contains something
            else:
                # Append the negative candidate to the valid negatives df
                valid_negatives = pd.concat([valid_negatives, negative_candidate], ignore_index=True)

                # If there are enough valid negatives (as many as the frequency of the miRNA in the miRNA family block)
                if len(valid_negatives) >= mirna_frequency:

                    # Get block rows for which column noncodingRNA == mirna and save to file (positives)
                    block_mirna = block[block['noncodingRNA'] == mirna].copy()
                    block_mirna.to_csv(output_file, sep='\t', index=False, header=False, mode='a')
                    
                    # Slice the valid negatives df to the required frequency, process valid negatives, and save to file (negatives)
                    valid_negatives = valid_negatives.iloc[:mirna_frequency].copy()
                    valid_negatives = process_valid_negatives(valid_negatives, block.columns)
                    valid_negatives.to_csv(output_file, sep='\t', index=False, header=False, mode='a')                    
                    # Exit the loop to move on to the next unique miRNA in the block
                    break
                # If there are not enough valid negatives    
                else:
                    # Check if the end of the negative gene pool is reached (edge case)
                    if index == len(negative_gene_pool) - 1:
                        if not valid_negatives.empty: 
                            block_mirna = block[block['noncodingRNA'] == mirna].iloc[:len(valid_negatives)].copy()
                            block_mirna.to_csv(output_file, sep='\t', index=False, header=False, mode='a')

                            valid_negatives = process_valid_negatives(valid_negatives, block.columns)
                            valid_negatives.to_csv(output_file, sep='\t', index=False, header=False, mode='a')

                            # Print a message if the end of the negative gene pool is reached and there are still not enough valid negatives
                            print(f"Missing {mirna_frequency - len(valid_negatives)} negatives for miRNA: {mirna}. Positive examples downsampled to retain 1:1 class ratio.", flush=True)
                        else:
                            # Print a message if the end of the negative gene pool is reached and there are no valid negatives
                            print(f"No valid negatives for miRNA: {mirna}. Excluding it from positive and negative examples.", flush=True)

def main():
    # Record start time
    start = time.time()

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Generate negative examples.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file name, MUST BE SORTED by 'miRNA family!'")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    parser.add_argument('--interaction_type', choices=['nonseed', 'canonicalseed', 'noncanonicalseed'], required=True, help="Interactions type to use for generating negative examples")
    args = parser.parse_args()
    
    # Read the entire positive examples file
    positive_samples = pd.read_csv(args.ifile, sep='\t')
    
    # Write header to the output file
    positive_samples.head(0).to_csv(args.ofile, sep='\t', index=False, mode='w')    

    # Get the set of all cluster ids
    all_clusters = positive_samples['gene_cluster_ID'].unique().tolist()

    with open(args.ofile, 'a') as ofile:
        
        # Iterate over blocks of positive examples with the same miRNA family
        for block in yield_mirnafam_blocks(args.ifile):

            # Check if the miRNA family is unknown
            if block['noncodingRNA_fam'].iloc[0] == 'unknown':
            
                # Get unique miRNA sequences (or names, equivalent)
                unique_mirnas = block['noncodingRNA'].unique().tolist()

                # Process each block of unique miRNA sequences
                for mirna in unique_mirnas:
                    
                    # Get the block of rows for this miRNA sequence
                    sub_block = block[block['noncodingRNA'] == mirna].copy()

                    # Run rest of code for each sub_block
                    process_block(sub_block, positive_samples, all_clusters, args.ofile, args.interaction_type)
                
                    print(f"Processed miRNA sequence block: {sub_block['noncodingRNA'].iloc[0]}", flush=True)

            else:
                # Process the block normally if miRNA family not 'unknown'
                process_block(block, positive_samples, all_clusters, args.ofile, args.interaction_type)

                print(f"Processed miRNA family block: {block['noncodingRNA_fam'].iloc[0]}", flush=True)
    
    # Record end time
    end = time.time() 

    # Print the difference between start and end time in secs
    print(f"The time of execution is:",(end-start),"s", flush=True)

if __name__ == "__main__":
    main()
