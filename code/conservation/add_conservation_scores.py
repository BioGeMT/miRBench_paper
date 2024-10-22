import pandas as pd
import numpy as np
import pyBigWig
import argparse
import warnings 

def read_input(input_file):
    # Read input data from a file or stdin.
    try:
        if input_file:
            data = pd.read_csv(input_file, sep='\t', dtype={'chr': str})
        else:
            data = pd.read_csv(sys.stdin, sep='\t', low_memory=False)
    
    except Exception as e:
        warnings.warn(f"Error reading input data: {e}", category=UserWarning)
        sys.exit(1)
    
    return data

# Check if a value can be converted to a float i.e. is a number, and if a number is not NaN
def is_valid_number(s):
    try:
        num = float(s)
        if np.isnan(num):
            return False
        return True

    except ValueError:
        return False

def get_valid_conservation(bw_file, chrom, start, end, chrom_sizes, gene):
    # Convert chrom to string
    chrom = str(chrom)
    
    # Special handling for mitochondrial chromosome
    if chrom == 'MT':
        chrom = 'chrM'
    elif not chrom.startswith('chr'):
        chrom = f"chr{chrom}"
    
    # Check if chromosome exists in the BigWig file
    if chrom not in chrom_sizes:
        print(f"Chromosome {chrom} not found in BigWig file. Returning NaN.")
        return np.nan

    try:
        # Convert start and end to integer, and adjust start to 0-based
        start = int(float(start)) - 1 
        end = int(float(end))

        # Get chromosome size
        chrom_size = chrom_sizes[chrom]

        # Check if interval is valid
        if start < 0 or start >= chrom_size or end <= 0 or end > chrom_size or start >= end:
            print(f"Invalid interval for {chrom}:{start}-{end}. Returning NaN.")
            return np.nan 

        # Get conservation score values
        conservation_scores = bw_file.values(chrom, start, end)

        # Check if all values are valid numbers
        if not all(is_valid_number(score) for score in conservation_scores):
            print(f"Invalid conservation score values for {chrom}:{start}-{end}. Returning NaN.")
            return np.nan

        # Check if length of conservation scores matches gene sequence length
        conservation_scores_length = len(conservation_scores)
        target_length = len(gene)
        if conservation_scores_length != target_length:
            print(f"Conservation scores length does not match gene sequence length for {chrom}:{start}-{end}. Returning NaN.")
            return np.nan

        # Return conservation scores if all checks pass    
        return conservation_scores

    except Exception as e:
        print(f"Error processing {chrom}:{start}-{end}: {e}. Returning NaN.")
        return np.nan

def add_conservation(df, phyloP_path, phastCons_path, ofile):
    # Open BigWig files
    with pyBigWig.open(phyloP_path) as bw_phyloP, pyBigWig.open(phastCons_path) as bw_phastCons:

        # Get chromosome sizes for both files
        chrom_sizes_phyloP = bw_phyloP.chroms()
        chrom_sizes_phastCons = bw_phastCons.chroms()

        with open(ofile, 'a') as ofile:                
            # Add conservation scores
            df['gene_phyloP'] = df.apply(lambda row: get_valid_conservation(bw_phyloP, row['chr'], row['start'], row['end'], chrom_sizes_phyloP, row['gene']), axis=1)
            df['gene_phastCons'] = df.apply(lambda row: get_valid_conservation(bw_phastCons, row['chr'], row['start'], row['end'], chrom_sizes_phastCons, row['gene']), axis=1)
                    
            df.to_csv(ofile, sep='\t', index=False, header=True, mode='a')

        return df

def main():
    parser = argparse.ArgumentParser(description="Add conservation scores.")
    parser.add_argument('--ifile', help="Input file (default: STDIN)")
    parser.add_argument('--ofile', help="Output file (default: STDOUT)")
    parser.add_argument('--phyloP_path', help="Path to phyloP BigWig file")
    parser.add_argument('--phastCons_path', help="Path to phastCons BigWig file")

    # Parse arguments
    args = parser.parse_args()

    # Read input data into a df
    df = read_input(args.ifile)

    # Add conservation scores
    add_conservation(df, args.phyloP_path, args.phastCons_path, args.ofile)

if __name__ == "__main__":
    main()