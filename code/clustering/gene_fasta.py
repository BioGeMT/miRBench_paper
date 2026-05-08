import pandas as pd
import argparse

def convert_tsv_to_fasta(input_file, output_file, lookup_file):
   
    data = pd.read_csv(input_file, sep='\t')
    unique_genes = data[['gene']].drop_duplicates().reset_index(drop=True)
    unique_genes['gene_id'] = range(1, len(unique_genes) + 1)
    unique_genes[['gene_id', 'gene']].to_csv(lookup_file, sep='\t', index=False)
    
    with open(output_file, 'w') as fasta_file:
        for index, row in unique_genes.iterrows():
            sequence = row['gene']
            fasta_file.write(f">{row['gene_id']}\n{sequence}\n")
    
    print(f"FASTA file created: {output_file}")
    print(f"Gene ID lookup file created: {lookup_file}")

def main():

    parser = argparse.ArgumentParser(description='Convert TSV file with gene sequences to FASTA format')
    parser.add_argument('--input', required=True, help='Input TSV file path')
    parser.add_argument('--output', required=True, help='Output FASTA file path')
    parser.add_argument('--lookup', required=True, help='Output TSV path for the gene ID lookup table')

    args = parser.parse_args()
    
    convert_tsv_to_fasta(args.input, args.output, args.lookup)

if __name__ == "__main__":
    main()
