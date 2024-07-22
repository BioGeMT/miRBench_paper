import argparse
import pandas as pd
import random
from Levenshtein import distance as levenshtein_distance

def generate_negative_samples(positive_samples, num_negatives, min_edit_distance):
    negative_samples = []

    # Create a list of tuples of unique mirna sequences and their corresponding families
    unique_smallRNAs_fams = list(positive_samples[['seq.m', 'noncodingRNA_fam']].drop_duplicates().values.tolist())
    gene_binding_dict = {}
    negative_blacklist = []
    unsuccessful = 0

    # Create a dictionary of lists of positive mirnas (values) binding to a gene (key)
    for _, row in positive_samples.iterrows():
        if str(row['seq.g']) in gene_binding_dict:
            if str(row['seq.m']) not in gene_binding_dict[str(row['seq.g'])]:
                gene_binding_dict[str(row['seq.g'])].append(str(row['seq.m']))
        else:
            gene_binding_dict[str(row['seq.g'])] = [str(row['seq.m'])]

    # Iterate over each positive sample to generate negatives
    for _, pos_row in positive_samples.iterrows():
        gene = pos_row['seq.g']
        n = num_negatives + unsuccessful
        for j in range(1, n + 1):
            distance = 0
            tries = 0
            while distance < min_edit_distance and tries < 200:
                random_mirna, random_fam = map(str, random.choice(unique_smallRNAs_fams))  
                distance = levenshtein_distance(random_mirna, gene_binding_dict[gene][0]) 
                for i in range(1, len(gene_binding_dict[gene])): 
                    d = levenshtein_distance(random_mirna, gene_binding_dict[gene][i]) 
                    distance = min(d, distance)
                if random_mirna+gene in negative_blacklist:
                    distance = 0
                tries += 1
            if tries == 200:
                unsuccessful = n - j
                break

            negative_sample = pos_row.copy()
            negative_sample['label'] = 0 
            negative_sample['seq.m'] = random_mirna
            negative_sample['noncodingRNA_fam'] = random_fam
            negative_samples.append(negative_sample)

            negative_blacklist.append(random_mirna+gene)
            unsuccessful = 0

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