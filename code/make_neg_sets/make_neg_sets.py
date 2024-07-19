import argparse
import random
import pandas as pd
from Levenshtein import distance as levenshtein_distance

def generate_negative_samples(positive_samples, num_negatives, min_edit_distance):
    negative_samples = []
    
    # Collect all unique seq.m and their corresponding noncodingRNA_fam values
    unique_seq_m_fams = positive_samples[['seq.m', 'noncodingRNA_fam']].drop_duplicates().values.tolist()
    
    gene_mirna_dict = {}
    for _, pos_row in positive_samples.iterrows():
        seq_g = pos_row['seq.g']
        seq_m = pos_row['seq.m']
        if seq_g not in gene_mirna_dict:
            gene_mirna_dict[seq_g] = set()
        gene_mirna_dict[seq_g].add(seq_m)

    negative_blacklist = set()
    unsuccessful = 0

    # Iterate over each positive sample to generate negatives
    for _, pos_row in positive_samples.iterrows():
        seq_g = pos_row['seq.g']
        
        # Calculate the number of negatives to generate for this positive sample
        # This includes making up for previously unsuccessful attempts
        num_negatives_for_sample = num_negatives + unsuccessful
        
        for i in range(num_negatives_for_sample):
            found_negative = False
            for attempt in range(200):  # Try up to 200 times to find a suitable negative
                random_mirna, random_fam = random.choice(unique_seq_m_fams)
                if (seq_g, random_mirna) in negative_blacklist:
                    continue
                
                distance = min(levenshtein_distance(random_mirna, pos_mirna) for pos_mirna in gene_mirna_dict[seq_g])
                
                if distance >= min_edit_distance:
                    negative_sample = pos_row.copy()
                    negative_sample['label'] = 0
                    negative_sample['seq.m'] = random_mirna
                    negative_sample['noncodingRNA_fam'] = random_fam  # Set the noncodingRNA_fam to the new seq.m's family
                    negative_samples.append(negative_sample)
                    negative_blacklist.add((seq_g, random_mirna))
                    found_negative = True
                    unsuccessful = 0
                    break
            
            if not found_negative:
                unsuccessful = num_negatives_for_sample - i - 1
                break

    return pd.DataFrame(negative_samples)

def main():
    parser = argparse.ArgumentParser(description="Generate negative samples with specific edit distance.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file name")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    parser.add_argument('--neg_ratio', type=int, default=100, help="Number of negative samples to generate per positive sample")
    parser.add_argument('--min_edit_distance', type=int, default=3, help="Minimum edit distance for negative samples")
    args = parser.parse_args()

    # Set fixed random seed
    random.seed(42)

    positive_samples = pd.read_csv(args.ifile, sep='\t')
    
    negative_samples = generate_negative_samples(positive_samples, args.neg_ratio, args.min_edit_distance)
    
    combined_data = pd.concat([positive_samples, negative_samples], ignore_index=True)
    
    combined_data.to_csv(args.ofile, sep='\t', index=False)

if __name__ == "__main__":
    main()

