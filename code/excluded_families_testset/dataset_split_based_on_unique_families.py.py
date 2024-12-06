import pandas as pd
import argparse

def filter_dataset(input_dataset, families_file, excluded_output, leftout_output):
    # Read the original dataset
    df = pd.read_csv(input_dataset, sep='\t')
    
    # Read the families file
    families_df = pd.read_csv(families_file, sep='\t')
    unique_families = set(families_df['noncodingRNA_fam'])
    
    # Split the dataset into two parts
    excluded_dataset = df[df['noncodingRNA_fam'].isin(unique_families)]
    remaining_dataset = df[~df['noncodingRNA_fam'].isin(unique_families)]
    
    # Save both datasets
    excluded_dataset.to_csv(excluded_output, sep='\t', index=False)
    remaining_dataset.to_csv(leftout_output, sep='\t', index=False)
    
    # Print some statistics
    print(f"Original dataset: {len(df)} rows")
    print(f"Excluded dataset: {len(excluded_dataset)} rows")
    print(f"Remaining dataset: {len(remaining_dataset)} rows")

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Original dataset TSV file')
parser.add_argument('--families', required=True, help='File with unique families')
parser.add_argument('--excluded_dataset', required=True, help='Output file for families that match input families file')
parser.add_argument('--remaining_dataset', required=True, help='Output file for families not in input families file')

args = parser.parse_args()
filter_dataset(args.input, args.families, args.excluded_dataset, args.remaining_dataset)
