import pandas as pd
import random
from Levenshtein import distance as levenshtein_distance
import argparse

def precompute_allowed_mirnas(positive_file_path, min_allowed_distance):
    '''
    precompute a dictionary of allowed mirnas for each mirna in the dataset.
    this speeds up the negative sample generation process.
    
     positive_file_path: path to the input file
     min_allowed_distance: minimum levenshtein distance for allowed mirnas
     dictionary of allowed mirnas for each mirna
    '''
    # read only necessary columns 
    df = pd.read_csv(positive_file_path, usecols=['noncodingRNA_seq', 'noncodingRNA_fam'])
    unique_mirnas = df.drop_duplicates().values.tolist()
    
    allowed_mirnas = {}
    for mirna1, fam1 in unique_mirnas:
        # store list of other mirnas (and their families) that are at least min_allowed_distance away
        allowed_mirnas[mirna1] = [(mirna2, fam2) for mirna2, fam2 in unique_mirnas 
                                  if mirna1 != mirna2 and levenshtein_distance(mirna1, mirna2) >= min_allowed_distance]
    return allowed_mirnas

def read_input_by_gene_batch(positive_file_handle):
    '''
    reads the input file in batches, grouped by gene.

    
     positive_file_handle: file handle for the input file
    yield: tuple of (header, block) where block is a list of rows for a single gene
    '''
    block = []
    gene = None
    header = next(positive_file_handle).strip().split('\t')
    
    for line in positive_file_handle:
        fields = line.strip().split('\t')
        current_gene = fields[header.index('seq.g')]
        
        if gene is None:
            gene = current_gene
        
        if current_gene == gene:
            block.append(fields)
        else:
            yield header, block
            block = [fields]
            gene = current_gene
    
    if block:
        yield header, block

def generate_negatives_gene(block, header, allowed_mirnas, num_negatives, min_required_edit_distance):
    '''
    generate negative samples for a single gene block.
    
     block: list of rows for a single gene
     header: list of column names
     allowed_mirnas: dictionary of allowed mirnas for each mirna
     num_negatives: number of negative samples to generate per positive sample
     min_required_edit_distance: minimum required levenshtein distance for negative samples
    return: list of negative samples
    '''
    negative_samples = []
    negative_blacklist = set()
    gene = block[0][header.index('seq.g')]
    mirna_index = header.index('noncodingRNA_seq')
    fam_index = header.index('noncodingRNA_fam')
    positive_mirnas = set(row[mirna_index] for row in block)
    
    for pos_row in block:
        current_pos_mirna = pos_row[mirna_index]
        # shuffle allowed mirnas to ensure randomness in selection
        shuffled_allowed_mirnas = random.sample(allowed_mirnas[current_pos_mirna], len(allowed_mirnas[current_pos_mirna]))
        
        negatives_found = 0
        for random_mirna, random_fam in shuffled_allowed_mirnas:
            if random_mirna not in negative_blacklist and random_mirna not in positive_mirnas:
                # double-check edit distance against all positive mirnas for the gene
                if all(levenshtein_distance(random_mirna, other_pos_mirna) >= min_required_edit_distance 
                       for other_pos_mirna in positive_mirnas if other_pos_mirna != current_pos_mirna):
                    neg_row = pos_row.copy()
                    neg_row[mirna_index] = random_mirna
                    neg_row[fam_index] = random_fam
                    neg_row[header.index('label')] = '0'
                    negative_samples.append(neg_row)
                    negative_blacklist.add(random_mirna)
                    negatives_found += 1
                    
                    if negatives_found == num_negatives:
                        break
        
        if negatives_found < num_negatives:
            print(f"warning: could only generate {negatives_found} out of {num_negatives} negatives for gene {gene} and mirna {current_pos_mirna}")
    
    return negative_samples

def main(positive_file_path, output_file_path, num_negatives, min_required_edit_distance):
    '''
    main function to generate negative samples for mirna-gene interactions.
    
     positive_file_path: path to the input file
     output_file_path: path to the output file
     num_negatives: number of negative samples to generate per positive sample
     min_required_edit_distance: minimum required levenshtein distance for negative samples
    '''
    # precompute allowed mirnas to speed up the process
    allowed_mirnas = precompute_allowed_mirnas(positive_file_path, min_required_edit_distance)
    
    with open(positive_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
        header_written = False
        # process input file in batches, grouped by gene
        for header, block in read_input_by_gene_batch(infile):
            if not header_written:
                outfile.write('\t'.join(header) + '\n')
                header_written = True
            
            # generate and write negative samples for each gene block
            negatives = generate_negatives_gene(block, header, allowed_mirnas, num_negatives, min_required_edit_distance)
            for neg_row in negatives:
                outfile.write('\t'.join(neg_row) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="generate negative samples for mirna-gene interactions.")
    parser.add_argument('--input', required=True, help="input file path")
    parser.add_argument('--output', required=True, help="output file path")
    parser.add_argument('--neg_ratio', type=int, default=1, help="number of negative samples to generate per positive sample")
    parser.add_argument('--min_edit_distance', type=int, default=3, help="minimum required edit distance for negative samples")
    args = parser.parse_args()

    main(args.input, args.output, args.neg_ratio, args.min_edit_distance)

     '''
it precomputes a list of allowed mirnas for each mirna based on a minimum edit distance. this speeds up the negative sample generation process.
it reads the input file in batches, grouped by gene. this allows processing of large files that might not fit in memory all at once.
for each gene, it generates negative samples:
it shuffles the list of allowed mirnas for each positive mirna.
it selects mirnas that are not in the positive set or already used as negatives.
it double-checks the edit distance against all positive mirnas for the gene to ensure diversity.
it creates negative samples up to the specified number or until it runs out of options.
it writes the generated negative samples to the output file.
if it can't generate the requested number of negatives for a gene-mirna pair, it prints a warning

details:

precompute_allowed_mirnas function:

this function prepares a dictionary of allowed mirnas for each mirna in the dataset.
it reads only the necessary columns (mirna sequence and family) from the input file to save memory.
it computes the levenshtein distance between each pair of mirnas.
for each mirna, it stores a list of other mirnas (and their families) that are at least min_allowed_distance away.
this precomputation step saves time during the negative sample generation process by reducing the search space.

read_input_by_gene_batch function:

this function implements a generator that reads the input file in batches, grouped by gene.
it reads the file line by line, grouping rows with the same gene together.
it yields each group (block) along with the header information.
this approach allows processing of very large files that might not fit in memory all at once.
it assumes the input file is sorted by gene. if not, a preprocessing step to sort the file would be necessary.

generate_negatives_gene function:

this is the core function that generates negative samples for each gene block.
it creates sets for positive mirnas and a blacklist for used negative mirnas to enable efficient lookups.
for each positive mirna in the block:
a. it shuffles the list of allowed mirnas to ensure randomness in selection.
b. it iterates through the shuffled list, looking for suitable negative samples.
c. it checks if the candidate mirna is not in the positive set or already used as a negative.
d. it performs an additional edit distance check against all positive mirnas for the gene.
e. if a suitable mirna is found, it creates a negative sample by copying the positive row and changing relevant fields.
f. it continues until it finds the requested number of negatives or exhausts the list of allowed mirnas.
if it can't generate enough negatives, it prints a warning message

'''

