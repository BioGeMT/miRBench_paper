import pandas as pd
import pyBigWig
import argparse
import warnings 

# Function to yield blocks of rows with the same gene from a sorted DataFrame
def yield_gene_blocks_from_df(df):
    current_block = []
    current_gene = None

    # Iterate over DataFrame rows
    for _, row in df.iterrows():
        gene_value = row['gene']
        
        if gene_value != current_gene:
            # If we have collected rows for the previous gene, yield them
            if current_block:
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

def add_inconsistent_column(df):

    df_inconsistent = pd.DataFrame(columns=['gene', 'noncodingRNA', 'noncodingRNA_fam', 'feature', 'test', 'label', 'chr', 'start', 'end', 'strand'])

    # Add 'inconsistent' column to DataFrame based on comparison of values within each gene block
    for gene_block in yield_gene_blocks_from_df(df):
        # Get the first row's values for comparison across the gene block
        chr_val = gene_block['chr'].iloc[0]
        start_val = gene_block['start'].iloc[0]
        end_val = gene_block['end'].iloc[0]
        strand_val = gene_block['strand'].iloc[0]

        # Add 'inconsistent' column based on comparison with the first row
        gene_block['inconsistent'] = gene_block.apply(
            lambda row: 1 if (row['chr'] != chr_val or 
                            row['start'] != start_val or 
                            row['end'] != end_val or 
                            row['strand'] != strand_val) else 0,
            axis=1
        )

        df_inconsistent = pd.concat([df_inconsistent, gene_block], ignore_index=True)

    return df_inconsistent

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

def add_conservation(df, phyloP_path, phastCons_path):
    # Open the BigWig files for phyloP and phastCons
    bw_phyloP = pyBigWig.open(phyloP_path)
    bw_phastCons = pyBigWig.open(phastCons_path)

    # Get chromosome sizes for both files
    chrom_sizes_phyloP = bw_phyloP.chroms()
    chrom_sizes_phastCons = bw_phastCons.chroms()

    # # Print unique chromosome names in input data
    # print("Unique chromosome names in input data:", df['chr'].unique())

    # Process each row for target sequence conservation
    df['target_phyloP'] = df.apply(lambda row: get_conservation(bw_phyloP, row['chr'], row['start'], row['end'], chrom_sizes_phyloP) if row['inconsistent'] == 1 else '-', axis=1)
    df['target_phastCons'] = df.apply(lambda row: get_conservation(bw_phastCons, row['chr'], row['start'], row['end'], chrom_sizes_phastCons) if row['inconsistent'] == 1 else '-', axis=1)

    # Close the BigWig files
    bw_phyloP.close()
    bw_phastCons.close()

    # # Check for empty results for phyloP scores
    # empty_results_phyloP = df[df['target_phyloP'].apply(lambda x: all(pd.isna(x)))]
    # if not empty_results_phyloP.empty:
    #     print(f"Number of rows with all NaN phyloP scores: {len(empty_results_phyloP)}")
    #     print("Sample of rows with all NaN phyloP scores:")
    #     print(empty_results_phyloP.head())
    #     print("Unique chromosomes in rows with NaN phyloP scores:", empty_results_phyloP['chr'].unique())

    # # Check for empty results for phastCons scores
    # empty_results_phastCons = df[df['target_phastCons'].apply(lambda x: all(pd.isna(x)))]
    # if not empty_results_phastCons.empty:
    #     print(f"Number of rows with all NaN phastCons scores: {len(empty_results_phastCons)}")
    #     print("Sample of rows with all NaN phastCons scores:")
    #     print(empty_results_phastCons.head())
    #     print("Unique chromosomes in rows with NaN phastCons scores:", empty_results_phastCons['chr'].unique())
    
    df = df.drop(columns=['inconsistent'])

    return df

# read input data from a file or stdin.
def read_input(input_file):

    try:
        if input_file:
            data = pd.read_csv(input_file, sep='\t', dtype={'chr': str})
        else:
            data = pd.read_csv(sys.stdin, sep='\t', low_memory=False)
    except Exception as e:
        warnings.warn(f"Error reading input data: {e}", category=UserWarning)
        sys.exit(1)
    
    return data

# write output data to a file or stdout.
def write_output(data, output_file):

    try:
        if output_file:
            data.to_csv(output_file, sep='\t', index=False)
        else:
            data.to_csv(sys.stdout, sep='\t', index=False)
    except Exception as e:
        warnings.warn(f"Error writing to output stream: {e}", category=UserWarning)
        sys.exit(1)

# main function to handle argument parsing and calling the processing functions.
def main():

    parser = argparse.ArgumentParser(description="Clean file after conservation scores have been added.")
    parser.add_argument('--ifile', help="Input file (default: STDIN)")
    parser.add_argument('--ofile', help="Output file (default: STDOUT)")
    parser.add_argument('--phyloP_path', help="Path to phyloP BigWig file")
    parser.add_argument('--phastCons_path', help="Path to phastCons BigWig file")

    args = parser.parse_args()

    # read input data.
    df = read_input(args.ifile)

    # add inconsistent column.
    df = add_inconsistent_column(df)

    # process data.
    df = add_conservation(df, args.phyloP_path, args.phastCons_path)

    # write output data.
    write_output(file_with_scores, args.ofile)

if __name__ == "__main__":
    main()