import pandas as pd
import pyBigWig
import argparse
import warnings 

# Function to yield blocks of rows with the same gene from a DataFrame (alphabetically) sorted by gene sequence
def yield_gene_blocks_from_df(df):
    current_block = []
    current_gene = None

    # Iterate over DataFrame rows
    for _, row in df.iterrows():
        gene_value = row['gene']
        
        if gene_value != current_gene:
            # Check if current_block contains rows for the previous gene and if so yield them
            if current_block: # returns True if not empty
                yield pd.DataFrame(current_block)
            # Start a new block for the new gene
            current_block = [row]
            current_gene = gene_value
        else:
            # If it's the same gene, keep adding to the current block
            current_block.append(row)

    # Yield the last block if it's not empty
    if current_block:
        yield pd.DataFrame(current_block)

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
    # Open the BigWig files for phyloP and phastCons
    bw_phyloP = pyBigWig.open(phyloP_path)
    bw_phastCons = pyBigWig.open(phastCons_path)

    # Get chromosome sizes for both files
    chrom_sizes_phyloP = bw_phyloP.chroms()
    chrom_sizes_phastCons = bw_phastCons.chroms()

    header = True
    columns = ['chr', 'start', 'end', 'strand']

    with open(ofile, 'a') as ofile:
        # Iterate over gene blocks
        for gene_block in yield_gene_blocks_from_df(df):
            
            # Add conservation scores to the gene block
            gene_block['target_phyloP'] = gene_block.apply(lambda row: get_conservation(bw_phyloP, row['chr'], row['start'], row['end'], chrom_sizes_phyloP), axis=1)
            gene_block['target_phastCons'] = gene_block.apply(lambda row: get_conservation(bw_phastCons, row['chr'], row['start'], row['end'], chrom_sizes_phastCons), axis=1)

            if header:
                gene_block.head(0).to_csv(ofile, sep='\t', index=False, mode='w')
                header = False
                
            gene_block.to_csv(ofile, sep='\t', index=False, header=False, mode='a')

    # Close the BigWig files
    bw_phyloP.close()
    bw_phastCons.close()
    

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