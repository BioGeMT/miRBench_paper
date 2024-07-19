import argparse
import random
import pandas as pd
from Levenshtein import distance as levenshtein_distance

def generate_negative_samples(positive_samples, num_negatives, min_edit_distance):
    negative_samples = []
    unique_seq_m_fams = positive_samples[['seq.m', 'noncodingRNA_fam']].drop_duplicates().values.tolist()
    
    gene_mirna_dict = {}
    for _, pos_row in positive_samples.iterrows():
        seq_g = pos_row['seq.g']
        seq_m = pos_row['seq.m']
        if seq_g not in gene_mirna_dict:
            gene_mirna_dict[seq_g] = set()
        gene_mirna_dict[seq_g].add(seq_m)

    generated_pairs = set()
    unsuccessful = 0

    for _, pos_row in positive_samples.iterrows():
        seq_g = pos_row['seq.g']
        n = num_negatives + unsuccessful
        
        for _ in range(n):
            potential_negatives = [
                (seq_m, noncodingRNA_fam) for seq_m, noncodingRNA_fam in unique_seq_m_fams
                if seq_m not in gene_mirna_dict[seq_g] and
                all(levenshtein_distance(seq_m, pos_seq_m) >= min_edit_distance 
                    for pos_seq_m in gene_mirna_dict[seq_g]) and
                (seq_g, seq_m) not in generated_pairs
            ]

            if not potential_negatives:
                unsuccessful = n - _ - 1
                break

            seq_m, noncodingRNA_fam = random.choice(potential_negatives)
            negative_sample = pos_row.copy()
            negative_sample['label'] = 0
            negative_sample['seq.m'] = seq_m
            negative_sample['noncodingRNA_fam'] = noncodingRNA_fam
            negative_samples.append(negative_sample)
            generated_pairs.add((seq_g, seq_m))
            unsuccessful = 0

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
