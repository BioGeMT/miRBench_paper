import pandas as pd
import argparse
import os
from seed_utils import find_seed_match

def process_dataset(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df['seed_type'] = df.apply(lambda row: find_seed_match(row['seq.g'], row['seq.m']), axis=1)
    
    # Count seed types for each miRNA family
    seed_counts = df.groupby('miRNA_fam')['seed_type'].value_counts().unstack(fill_value=0)
    
    # Ensure all seed type columns exist
    for seed_type in ['Seed6mer', 'Seed7mer', 'Seed8mer']:
        if seed_type not in seed_counts.columns:
            seed_counts[seed_type] = 0
    
    # Calculate non-canonical seeds
    seed_counts['SeedNoncanonical'] = seed_counts.sum(axis=1) - seed_counts[['Seed6mer', 'Seed7mer', 'Seed8mer']].sum(axis=1)
    
    # Calculate total entries for the dataset
    total_entries = seed_counts.sum().sum()
    
    # Calculate percentages
    seed_percentages = seed_counts / total_entries * 100
    
    # Reorder columns
    counts = seed_counts[['Seed6mer', 'Seed7mer', 'Seed8mer', 'SeedNoncanonical']]
    percentages = seed_percentages[['Seed6mer', 'Seed7mer', 'Seed8mer', 'SeedNoncanonical']]
    
    return counts, percentages

def main(input_files, output_file):
    all_counts = []
    all_percentages = []
    for file_path in input_files:
        counts, percentages = process_dataset(file_path)
        dataset_name = os.path.splitext(os.path.basename(file_path))[0]
        
        counts.columns = [f'{dataset_name}_{col[4:].lower()}_count' for col in counts.columns]
        percentages.columns = [f'{dataset_name}_{col[4:].lower()}_percent' for col in percentages.columns]
        
        all_counts.append(counts)
        all_percentages.append(percentages)
    
    # Combine all datasets
    result_counts = pd.concat(all_counts, axis=1)
    result_percentages = pd.concat(all_percentages, axis=1)
    
    # Combine counts and percentages
    result = pd.concat([result_counts, result_percentages], axis=1)
    
    # Calculate total count for sorting
    result['total_count'] = result[[col for col in result.columns if col.endswith('_count')]].sum(axis=1)
    
    # Sort by total count in descending order
    result = result.sort_values('total_count', ascending=False)
    
    # Remove the total_count column
    result = result.drop('total_count', axis=1)
    
    # Reset index to get miRNA family as a column
    result = result.reset_index()
    
    # Handle NaN values and convert to appropriate types
    for col in result.columns:
        if col.endswith('_percent'):
            result[col] = result[col].fillna(0).round(2)
        elif col.endswith('_count'):
            result[col] = result[col].fillna(0).astype(int)
    
    # Save the result to TSV file
    result.to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to {output_file}")
    
    # Print summary information
    print("\nDataframe shape:", result.shape)
    print("\nColumn names:", result.columns.tolist())
    print("\nFirst few rows of the dataframe:")
    print(result.head())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process miRNA datasets and create a combined seed type analysis with counts and percentages.')
    parser.add_argument('input_files', nargs='+', help='Input TSV files for datasets')
    parser.add_argument('-o', '--output', default='mirna_seed_analysis_counts_percentages.tsv', help='Output TSV file (default: mirna_seed_analysis_counts_percentages.tsv)')
    args = parser.parse_args()
    
    main(args.input_files, args.output)