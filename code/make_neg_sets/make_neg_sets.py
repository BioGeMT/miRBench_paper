import argparse
import random
import pandas as pd
from Levenshtein import distance as levenshtein_distance

def generate_negative_samples(positive_samples, num_negatives, min_edit_distance):
    negative_samples = []
    unique_seq_m = positive_samples['seq.m'].unique()
    
    gene_mirna_dict = {}
    for _, pos_row in positive_samples.iterrows():
        seq_g = pos_row['seq.g']
        seq_m = pos_row['seq.m']
        if seq_g not in gene_mirna_dict:
            gene_mirna_dict[seq_g] = []
        gene_mirna_dict[seq_g].append(seq_m)

    negative_blacklist = set()
    unsuccessful = 0

    for seq_g, positive_mirnas in gene_mirna_dict.items():
        n = len(positive_mirnas) * num_negatives + unsuccessful
        for _ in range(n):
            found_negative = False
            for _ in range(200):  # Try up to 200 times
                random_mirna = random.choice(unique_seq_m)
                if (seq_g, random_mirna) in negative_blacklist:
                    continue
                
                distance = min(levenshtein_distance(random_mirna, pos_mirna) for pos_mirna in positive_mirnas)
                
                if distance >= min_edit_distance:
                    negative_sample = positive_samples[positive_samples['seq.g'] == seq_g].iloc[0].copy()
                    negative_sample['label'] = 0
                    negative_sample['seq.m'] = random_mirna
                    negative_samples.append(negative_sample)
                    negative_blacklist.add((seq_g, random_mirna))
                    found_negative = True
                    unsuccessful = 0
                    break
            
            if not found_negative:
                unsuccessful = n - _ - 1
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

