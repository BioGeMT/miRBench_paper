import pandas as pd
import argparse

def filter_dataset(input_dataset, families_file, excluded_output, leftout_output):
    # Read the original dataset
    df = pd.read_csv(input_dataset, sep='\t')
    
    # Read the families file (skipping the 'count' column)
    families_df = pd.read_csv(families_file, sep='\t')
    unique_families = set(families_df['noncodingRNA_fam'])
    
    # Split the dataset into two parts
    excluded_df = df[df['noncodingRNA_fam'].isin(unique_families)]
    leftout_df = df[~df['noncodingRNA_fam'].isin(unique_families)]
    
    # Save both datasets
    excluded_df.to_csv(excluded_output, sep='\t', index=False)
    leftout_df.to_csv(leftout_output, sep='\t', index=False)
    
    # Print some statistics
    print(f"Original dataset: {len(df)} rows")
    print(f"Excluded dataset: {len(excluded_df)} rows")
    print(f"Leftout dataset: {len(leftout_df)} rows")

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Original dataset TSV file')
parser.add_argument('--families', required=True, help='File with unique families')
parser.add_argument('--excluded', required=True, help='Output file for families that match input families file')
parser.add_argument('--leftout', required=True, help='Output file for families not in input families file')

args = parser.parse_args()
filter_dataset(args.input, args.families, args.excluded, args.leftout)
