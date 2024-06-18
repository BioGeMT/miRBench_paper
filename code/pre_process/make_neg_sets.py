import argparse
import random
import pandas as pd
from Levenshtein import distance as levenshtein_distance

def generate_negative_samples(data, positive_samples, num_negatives, min_edit_distance):
    negative_samples = []

    # Collect all unique seq.m and their corresponding noncodingRNA_fam values
    unique_seq_m_fams = data[['seq.m', 'noncodingRNA_fam']].drop_duplicates()

    for _, pos_row in positive_samples.iterrows():
        seq_g = pos_row['seq.g']
        seq_m_original = pos_row['seq.m']

        # Filter out the seq.m values that are already used in the positive samples with the same seq.g
        potential_negatives = [
            (seq_m, noncodingRNA_fam) for seq_m, noncodingRNA_fam in unique_seq_m_fams.values
            if seq_m != seq_m_original and levenshtein_distance(seq_m_original, seq_m) >= min_edit_distance
        ]

        # If not enough potential negatives found, raise an error or continue with fewer negatives
        if len(potential_negatives) < num_negatives:
            print(f"Warning: Only {len(potential_negatives)} potential negative samples found for seq.g: {seq_g}")
            random.shuffle(potential_negatives)
        else:
            random.shuffle(potential_negatives)
            potential_negatives = potential_negatives[:num_negatives]

        for seq_m, noncodingRNA_fam in potential_negatives:
            negative_sample = pos_row.copy()
            negative_sample['label'] = 0
            negative_sample['seq.m'] = seq_m
            negative_sample['noncodingRNA_fam'] = noncodingRNA_fam  # Set the noncodingRNA_fam to the new seq.m's family
            negative_samples.append(negative_sample)

    return pd.DataFrame(negative_samples)

def main():
    parser = argparse.ArgumentParser(description="Generate negative samples with specific edit distance.")
    parser.add_argument('--ifile', type=str, required=True, help="Input file name")
    parser.add_argument('--ofile', type=str, required=True, help="Output file name")
    parser.add_argument('--neg_ratio', type=int, default=100, help="Number of negative samples to generate per positive sample")
    parser.add_argument('--min_edit_distance', type=int, default=3, help="Minimum edit distance for negative samples")

    args = parser.parse_args()

    # Load the input data
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
