import argparse
import pandas as pd
import random
from Levenshtein import distance as levenshtein_distance

def generate_negative_samples(positive_samples, num_negatives, min_required_edit_distance):
    """
    Generate negative samples based on positive samples for a miRNA-gene interaction dataset.

    This function creates negative samples by randomly selecting miRNAs that are not known
    to interact with specific genes, while ensuring a minimum edit distance from known 
    positive interactions.

    Parameters:
    positive_samples (pd.DataFrame): A DataFrame containing positive miRNA-gene interactions.
                                     Must include columns 'seq.g' (gene), 'seq.m' (miRNA),
                                     and 'noncodingRNA_fam' (miRNA family). Should also include
                                     a 'label' column with value 1 for positive samples.
    num_negatives (int): The number of negative samples to generate for each positive sample.
    min_required_edit_distance (int): The minimum edit distance required between a negative
                                      sample's miRNA and all positive miRNAs for the same gene.

    Returns:
    pd.DataFrame: A DataFrame containing the generated negative samples, with the same
                  structure as the input positive_samples, but with 'label' set to 0.

    Notes:
    - The function attempts to maintain a balanced dataset by generating 'num_negatives'
      negative samples for each positive sample.
    - If the function fails to generate all required negative samples for a gene, it will
      attempt to compensate in subsequent iterations.
    - A maximum of 200 attempts are made to find a suitable negative sample for each required
      negative interaction.
    - The function uses a blacklist to avoid generating duplicate negative samples.
    - The function uses the Levenshtein distance to ensure a minimum edit distance between 
        the potential miRNA for the negative sample and all the positive miRNAs for the same gene.
    - If all attempts to generate negative samples fail, a warning is printed.
    """

    negative_samples = []

    # Create a list of tuples of unique mirna sequences and their corresponding families
    unique_smallRNAs_fams = list(positive_samples[['seq.m', 'noncodingRNA_fam']].drop_duplicates().values.tolist())
    gene_binding_dict = {}
    negative_blacklist = []
    unsuccessful = 0

    # Create a dictionary of lists of unique positive mirnas (values) binding to a gene (key)
    for _, row in positive_samples.iterrows():
        if str(row['seq.g']) in gene_binding_dict:
            if str(row['seq.m']) not in gene_binding_dict[str(row['seq.g'])]:
                gene_binding_dict[str(row['seq.g'])].append(str(row['seq.m']))
        else:
            gene_binding_dict[str(row['seq.g'])] = [str(row['seq.m'])]

    # Iterate over each positive sample to generate negatives
    for i, pos_row in positive_samples.iterrows():
        gene = pos_row['seq.g']
        # The number of negative samples to generate for this positive sample is the sum of the number of negative samples to generate and the number of remaining negative samples to generate due to a failed attempt in the previous positive sample
        n = num_negatives + unsuccessful
        for j in range(1, n + 1): # Iterating over 1 to n+1 to easily account for unsuccessful attempts later
            min_edit_distance = 0
            tries = 0
            while min_edit_distance < min_required_edit_distance and tries < 200:
                random_mirna, random_fam = map(str, random.choice(unique_smallRNAs_fams))  
                # Check if the random mirna is not already binding to the gene
                if random_mirna in gene_binding_dict[gene]:
                    continue
                # Check if the random mirna is not already in the blacklist
                if random_mirna + gene in negative_blacklist:
                    continue
                # Get a minimum edit distance from all the positive mirnas of given gene
                min_edit_distance = min([levenshtein_distance(random_mirna, pos_mirna) for pos_mirna in gene_binding_dict[gene]])
                tries += 1
            if tries == 200:
                # Update `unsuccessful` to add the remaining number of negative samples to generate for this positive sample to the next positive sample to maintain the positive to negative ratio
                unsuccessful = n - j
                if i == len(positive_samples.index) - 1:
                    print(f"Warning: Failed to generate all negative samples. Missing {unsuccessful} negative samples.")
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
    parser.add_argument('--min_required_edit_distance', type=int, default=3, help="Minimum required edit distance for negative samples")
    args = parser.parse_args()

    # Set a fixed random seed for reproducibility
    random.seed(42)

    # Load positive samples from the input file
    positive_samples = pd.read_csv(args.ifile, sep='\t')
    
    # Generate negative samples
    negative_samples = generate_negative_samples(positive_samples, args.neg_ratio, args.min_required_edit_distance)
    
    # Combine positive and negative samples
    combined_data = pd.concat([positive_samples, negative_samples], ignore_index=True)
    
    # Save the combined dataset to the output file
    combined_data.to_csv(args.ofile, sep='\t', index=False)

if __name__ == "__main__":
    main()