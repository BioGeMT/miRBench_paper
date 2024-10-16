import pandas as pd
import pyBigWig
import argparse
import warnings 

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
        return [float('nan')] * (int(end) - int(start))

    try:
        start = int(float(start)) - 1  # Convert to 0-based
        end = int(float(end))
        chrom_size = chrom_sizes[chrom]
        if start >= chrom_size or end <= 0 or start >= end:
            print(f"Invalid interval for {chrom}:{start}-{end}")
            return [float('nan')] * (end - start)
        start = max(0, start)
        end = min(chrom_size, end)
        return bw_file.values(chrom, start, end)
    except Exception as e:
        print(f"Error processing {chrom}:{start}-{end}: {e}")
        return [float('nan')] * (end - start)

def add_conservation(df, phyloP_path, phastCons_path, ofile):
    # Open BigWig files
    with pyBigWig.open(phyloP_path) as bw_phyloP, pyBigWig.open(phastCons_path) as bw_phastCons:

        # Get chromosome sizes for both files
        chrom_sizes_phyloP = bw_phyloP.chroms()
        chrom_sizes_phastCons = bw_phastCons.chroms()

        header = True
        columns = ['chr', 'start', 'end', 'strand']

        with open(ofile, 'a') as ofile:                
            # Add conservation scores
            df['gene_phyloP'] = df.apply(lambda row: get_conservation(bw_phyloP, row['chr'], row['start'], row['end'], chrom_sizes_phyloP), axis=1)
            df['gene_phastCons'] = df.apply(lambda row: get_conservation(bw_phastCons, row['chr'], row['start'], row['end'], chrom_sizes_phastCons), axis=1)

            if header:
                df.head(0).to_csv(ofile, sep='\t', index=False, mode='w')
                header = False
                    
            df.to_csv(ofile, sep='\t', index=False, header=False, mode='a')
    

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