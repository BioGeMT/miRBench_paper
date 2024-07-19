import argparse
import random
import pandas as pd
from Levenshtein import distance as levenshtein_distance

def generate_negative_samples(positive_samples, num_negatives, min_edit_distance):
    negative_samples = []
    
    # Create a list of all miRNA sequences and their corresponding RNA families
    all_seq_m_fams = positive_samples[['seq.m', 'noncodingRNA_fam']].values.tolist()
    
    # Create a dictionary mapping each gene to its set of positive miRNAs
    gene_mirna_dict = {}
    for _, pos_row in positive_samples.iterrows():
        seq_g = pos_row['seq.g']
        seq_m = pos_row['seq.m']
        if seq_g not in gene_mirna_dict:
            gene_mirna_dict[seq_g] = set()
        gene_mirna_dict[seq_g].add(seq_m)

    # Set to keep track of gene-miRNA pairs already used as negative samples
    negative_blacklist = set()

    # Counter for unsuccessful attempts
    unsuccessful = 0

    # Iterate over each positive sample to generate negatives
    for _, pos_row in positive_samples.iterrows():
        seq_g = pos_row['seq.g']
        
        # Calculate the number of negatives to generate for this positive sample
        # This includes making up for previously unsuccessful attempts
        num_negatives_for_sample = num_negatives + unsuccessful
        
        negatives_generated = 0
        for _ in range(num_negatives_for_sample):
            found_negative = False
            # Try up to 200 times to find a suitable negative sample
            for attempt in range(200):
                # Randomly select a miRNA and its family from all entries
                random_mirna, random_fam = random.choice(all_seq_m_fams)
                # Check if this gene-miRNA pair has already been used as a negative sample
                if (seq_g, random_mirna) in negative_blacklist:
                    continue
                
                # Check if the selected miRNA is sufficiently different from all positive miRNAs for this gene
                if all(levenshtein_distance(random_mirna, pos_mirna) >= min_edit_distance 
                       for pos_mirna in gene_mirna_dict[seq_g]):
                    # Create a new negative sample
                    negative_sample = pos_row.copy()
                    negative_sample['label'] = 0
                    negative_sample['seq.m'] = random_mirna
                    negative_sample['noncodingRNA_fam'] = random_fam
                    negative_samples.append(negative_sample)
                    # Add this gene-miRNA pair to the blacklist
                    negative_blacklist.add((seq_g, random_mirna))
                    found_negative = True
                    negatives_generated += 1
                    break
            
            if not found_negative:
                break

        # Update the unsuccessful counter
        unsuccessful = num_negatives_for_sample - negatives_generated

        if unsuccessful > 0:
            print(f"Warning: Could not generate all {num_negatives_for_sample} negative samples for gene {seq_g}. Generated {negatives_generated}. Unsuccessful: {unsuccessful}")

    return pd.DataFrame(negative_samples)

def main():
    parser = argparse.ArgumentParser(description="Generate negative samples with specific edit distance.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file name")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    parser.add_argument('--neg_ratio', type=int, default=100, help="Number of negative samples to generate per positive sample")
    parser.add_argument('--min_edit_distance', type=int, default=3, help="Minimum edit distance for negative samples")
    args = parser.parse_args()

    # Set a fixed random seed for reproducibility
    random.seed(42)

    # Load positive samples from the input file
    positive_samples = pd.read_csv(args.ifile, sep='\t')
    
    # Generate negative samples
    negative_samples = generate_negative_samples(positive_samples, args.neg_ratio, args.min_edit_distance)
    
    # Combine positive and negative samples
    combined_data = pd.concat([positive_samples, negative_samples], ignore_index=True)
    
    # Save the combined dataset to the output file
    combined_data.to_csv(args.ofile, sep='\t', index=False)

if __name__ == "__main__":
    main()

