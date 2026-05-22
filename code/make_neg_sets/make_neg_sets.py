import argparse
import pandas as pd
import time
import hashlib

def yield_negative_sampling_sub_blocks(block):
    # Check if the miRNA family is unknown
    if block['noncodingRNA_fam'].iloc[0] == 'unknown':
    
        # Get unique miRNA sequences (or names, equivalent)
        unique_mirnas = block['noncodingRNA'].unique().tolist()

        # Process each block of unique miRNA sequences
        for mirna in unique_mirnas:
            
            # Get the block of rows for this miRNA sequence
            sub_block = block[block['noncodingRNA'] == mirna]

            yield sub_block

    else:
        # Process the block normally if not 'unknown'
        yield block

def get_negative_sampling_block_label(block):
    if block['noncodingRNA_fam'].iloc[0] == 'unknown':
        return block['noncodingRNA'].iloc[0]
    return block['noncodingRNA_fam'].iloc[0]

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

def process_block(block, positive_samples, all_clusters, output_file):

    # Set a fixed seed for reproducibility
    ## Get the first item from 'noncodingRNA_name'
    miRNA_name = block['noncodingRNA_name'].iloc[0]

    ## Generate a SHA-256 hash and get the hexadecimal string
    miRNA_hash_hex = hashlib.sha256(miRNA_name.encode()).hexdigest()

    ## Convert the hexadecimal hash to a decimal integer
    miRNA_hash_int = int(miRNA_hash_hex, 16)

    ## Reduce the size using modulo (e.g., within the range of a 32-bit unsigned integer) and set it as the seed
    seed = miRNA_hash_int % 4294967295

    # Get the set of cluster ids that share this miRNA family
    block_clusters = block['gene_cluster_ID'].unique().tolist()

    # Get the set of cluster ids that are allowed to be paired with this miRNA family
    mirfam_allowed_clusters = [cluster for cluster in all_clusters if cluster not in block_clusters]

    # Pool gene rows from allowed clusters
    negative_pool = positive_samples[positive_samples['gene_cluster_ID'].isin(mirfam_allowed_clusters)]

    # Shuffle the negative pool and drop duplicates based on ClusterID
    negative_pool = negative_pool.sample(frac=1, random_state=seed).drop_duplicates(subset=['gene_cluster_ID'], keep='first')

    # Get the number of negatives to be generated for this miRNA family block
    num_neg = block.shape[0]
    block_label = get_negative_sampling_block_label(block)

    if len(negative_pool) == 0:
        print(
            f"Warning: Skipping block {block_label} - {block.shape[0]} positives excluded (no eligible negative clusters available)",
            flush=True,
        )
        return

    if num_neg > len(negative_pool):
        excluded = num_neg - len(negative_pool)
        print(
            f"Warning: Not enough negative examples for block {block_label}. {excluded} positives excluded",
            flush=True,
        )
        block = block.sample(n=len(negative_pool), random_state=seed)
        num_neg = len(negative_pool)

    # Sample num_neg from mirfam_allowed_genes rows
    negative_genes = negative_pool.sample(n=num_neg, random_state=seed)

    # Start constructing the df rows for the negative examples
    columns = ['gene', 'feature', 'test', 'chr', 'start', 'end', 'strand', 'gene_cluster_ID', 'dominant_region', 'regions_present', 'read_start_in_sel_tx_1based', 'read_end_in_sel_tx_1based']
    negatives_df  = negative_genes[columns].copy()

    # Add the miRNA sequence, name and family columns from block to negatives_df by index
    negatives_df['noncodingRNA'] = block['noncodingRNA'].values
    negatives_df['noncodingRNA_name'] = block['noncodingRNA_name'].values
    negatives_df['noncodingRNA_fam'] = block['noncodingRNA_fam'].values

    # Add the label column to negatives_df
    negatives_df['label'] = 0

    # Add Nunique column set to 0 for negative examples
    negatives_df['Nunique'] = 0

    # Reorder columns in negatives_df to match the order in block
    negatives_df = negatives_df[block.columns]

    # Append positive examples for this block to the output file
    block.to_csv(output_file, sep='\t', index=False, header=False, mode='a')

    # Append negative examples for this block to the output file
    negatives_df.to_csv(output_file, sep='\t', index=False, header=False, mode='a')

def main():
    # Record start time
    start = time.time()

    parser = argparse.ArgumentParser(description="Generate negative examples.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file name, MUST BE SORTED by 'miRNA family!'")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    args = parser.parse_args()
    
    # Read the entire positive examples file
    positive_samples = pd.read_csv(args.ifile, sep='\t')
    
    # Write header to the output file
    positive_samples.head(0).to_csv(args.ofile, sep='\t', index=False, mode='w')    

    # Get the set of all cluster ids
    all_clusters = positive_samples['gene_cluster_ID'].unique().tolist()

    with open(args.ofile, 'a') as ofile:
        
        for block in yield_mirnafam_blocks(args.ifile):

            for sub_block in yield_negative_sampling_sub_blocks(block):

                # Run rest of code for each sub_block
                process_block(sub_block, positive_samples, all_clusters, args.ofile)
            
                if sub_block['noncodingRNA_fam'].iloc[0] == 'unknown':
                    print(f"Processed miRNA sequence block: {get_negative_sampling_block_label(sub_block)}", flush=True)
                else:
                    print(f"Processed miRNA family block: {get_negative_sampling_block_label(sub_block)}", flush=True)
    
    # Record end time
    end = time.time() 

    # Print the difference between start and end time in secs
    print(f"The time of execution is:",(end-start),"s", flush=True)

if __name__ == "__main__":
    main()
