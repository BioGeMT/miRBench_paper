import argparse
import random
import pandas as pd
from Levenshtein import distance as levenshtein_distance

def generate_negative_samples(data, positive_samples, num_negatives, min_edit_distance):
    negative_samples = []
    # Get unique seq.m and noncodingRNA_fam pairs from the entire dataset
    unique_seq_m_fams = data[['seq.m', 'noncodingRNA_fam']].drop_duplicates().values.tolist()
    
    # Create a dictionary to store all positive miRNAs for each gene
    gene_mirna_dict = {}
    for _, pos_row in positive_samples.iterrows():
        seq_g = pos_row['seq.g']
        seq_m = pos_row['seq.m']
        if seq_g not in gene_mirna_dict:
            gene_mirna_dict[seq_g] = set()
        gene_mirna_dict[seq_g].add(seq_m)

    # Set to keep track of generated negative pairs
    generated_pairs = set()

    # Iterate through each positive sample to generate negative samples
    for _, pos_row in positive_samples.iterrows():
        seq_g = pos_row['seq.g']
        
        # Find potential negative samples for this seq_g
        potential_negatives = [
            (seq_m, noncodingRNA_fam) for seq_m, noncodingRNA_fam in unique_seq_m_fams
            if seq_m not in gene_mirna_dict[seq_g] and  # Ensure seq_m is not a positive for this seq_g
            all(levenshtein_distance(seq_m, pos_seq_m) >= min_edit_distance 
                for pos_seq_m in gene_mirna_dict[seq_g]) and  # Check edit distance against all positives
            (seq_g, seq_m) not in generated_pairs  # Ensure this pair hasn't been generated before
        ]

        # Warn if not enough potential negatives found
        if len(potential_negatives) < num_negatives:
            print(f"Warning: Only {len(potential_negatives)} potential negative samples found for seq.g: {seq_g}")
        
        # Randomly select negative samples
        selected_negatives = random.sample(potential_negatives, min(num_negatives, len(potential_negatives)))
        
        # Create negative samples
        for seq_m, noncodingRNA_fam in selected_negatives:
            negative_sample = pos_row.copy()
            negative_sample['label'] = 0  # Set label to 0 for negative samples
            negative_sample['seq.m'] = seq_m
            negative_sample['noncodingRNA_fam'] = noncodingRNA_fam
            negative_samples.append(negative_sample)
            generated_pairs.add((seq_g, seq_m))  # Mark this pair as generated

    return pd.DataFrame(negative_samples)

def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description="Generate negative samples with specific edit distance.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file name")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    parser.add_argument('--neg_ratio', type=int, default=100, help="Number of negative samples to generate per positive sample")
    parser.add_argument('--min_edit_distance', type=int, default=3, help="Minimum edit distance for negative samples")
    args = parser.parse_args()

    # Read input data
    data = pd.read_csv(args.ifile, sep='\t')
    
    # Filter positive samples
    positive_samples = data[data['label'] == 1]
    
    # Generate negative samples
    negative_samples = generate_negative_samples(data, positive_samples, args.neg_ratio, args.min_edit_distance)
    
    # Combine positive and negative samples
    combined_data = pd.concat([data, negative_samples], ignore_index=True)
    
    # Save the result to a new file
    combined_data.to_csv(args.ofile, sep='\t', index=False)

if __name__ == "__main__":
    main()
