import pandas as pd
import argparse

def filter_dataset(input_dataset, families_file, output_path):
    # Read the original dataset
    df = pd.read_csv(input_dataset, sep='\t')
    
    # Read the families file (skipping the 'count' column)
    families_df = pd.read_csv(families_file, sep='\t')
    unique_families = set(families_df['noncodingRNA_fam'])
    
    # Filter the original dataset to keep only rows with the unique families
    filtered_df = df[df['noncodingRNA_fam'].isin(unique_families)]
    
    # Save the filtered dataset
    filtered_df.to_csv(output_path, sep='\t', index=False)

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Original dataset TSV file')
parser.add_argument('--families', required=True, help='File with unique families')
parser.add_argument('--output', required=True, help='Output filtered dataset')

args = parser.parse_args()
filter_dataset(args.input, args.families, args.output)