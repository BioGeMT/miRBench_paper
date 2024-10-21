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

# Check if a value can be converted to a float i.e. is a valid number
def is_valid_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def get_conservation(bw_file, chrom, start, end, chrom_sizes):
    # Convert to string
    chrom = str(chrom)
    
    # Special handling for mitochondrial chromosome
    if chrom == 'MT':
        chrom = 'chrM'
    elif not chrom.startswith('chr'):
        chrom = f"chr{chrom}"
    
    # Check if chromosome exists in the BigWig file
    if chrom not in chrom_sizes:
        print(f"Chromosome {chrom} not found in BigWig file")
        return np.nan

    try:
        start = int(float(start)) - 1  # Convert to 0-based
        end = int(float(end))
        chrom_size = chrom_sizes[chrom]
        # Check if interval is valid
        if start < 0 or start >= chrom_size or end <= 0 or end > chrom_size or start >= end:
            print(f"Invalid interval for {chrom}:{start}-{end}")
            return np.nan 
        # Get conservation score values
        conservation_scores = bw_file.values(chrom, start, end)
        # Check if all values are valid numbers
        if all(is_valid_number(score) for score in conservation_scores):
            return conservation_scores
        else:
            print(f"Invalid conservation score numbers for {chrom}:{start}-{end}")
            return np.nan
    except Exception as e:
        print(f"Error processing {chrom}:{start}-{end}: {e}")
        return np.nan

def check_conservation_length(row, column):
    try:
        conservation_scores_length = len(row[column])
        target_length = len(row['gene'])
        if conservation_scores_length != target_length:
            print(f"Conservation scores length does not match gene sequence length for {row['chr']}:{row['start']}-{row['end']}")
            return False
        else:
            return True
    except:
        print(f"Error checking conservation scores length for {row['chr']}:{row['start']}-{row['end']}")
        return False

def add_conservation(df, phyloP_path, phastCons_path, ofile):
    # Open BigWig files
    with pyBigWig.open(phyloP_path) as bw_phyloP, pyBigWig.open(phastCons_path) as bw_phastCons:

        # Get chromosome sizes for both files
        chrom_sizes_phyloP = bw_phyloP.chroms()
        chrom_sizes_phastCons = bw_phastCons.chroms()

        columns = ['chr', 'start', 'end', 'strand']

        with open(ofile, 'a') as ofile:                
            # Add conservation scores
            df['gene_phyloP'] = df.apply(lambda row: get_conservation(bw_phyloP, row['chr'], row['start'], row['end'], chrom_sizes_phyloP), axis=1)
            df['gene_phastCons'] = df.apply(lambda row: get_conservation(bw_phastCons, row['chr'], row['start'], row['end'], chrom_sizes_phastCons), axis=1)

            # Check if conservation scores have the same length as the gene sequence, otherwise set to NaN
            df['gene_phyloP'] = df.apply(lambda row: row['gene_phyloP'] if check_conservation_length(row, 'gene_phyloP') else np.nan, axis=1)
            df['gene_phastCons'] = df.apply(lambda row: row['gene_phastCons'] if check_conservation_length(row, 'gene_phastCons') else np.nan, axis=1)
                    
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