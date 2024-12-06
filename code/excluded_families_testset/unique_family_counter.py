import pandas as pd
import argparse

def read_and_filter_df(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df[(df['noncodingRNA_fam'] != 'unknown') & (df['noncodingRNA_fam'] != '0')]

def analyze_unique_families(unique_input_file, file2_path, file3_path, output_path):
    df1 = read_and_filter_df(unique_input_file)
    df2 = read_and_filter_df(file2_path)
    df3 = read_and_filter_df(file3_path)
    
    families1 = set(df1['noncodingRNA_fam'].unique())
    families2 = set(df2['noncodingRNA_fam'].unique())
    families3 = set(df3['noncodingRNA_fam'].unique())
    
    unique_to_1 = families1 - (families2 | families3)
    unique_counts = df1[df1['noncodingRNA_fam'].isin(unique_to_1)]['noncodingRNA_fam'].value_counts()
    
    unique_counts.to_csv(output_path, sep='\t')

    print(f"Unique families: {len(unique_to_1)}")
    print(f"Total occurrences: {unique_counts.sum()}")

parser = argparse.ArgumentParser()
parser.add_argument('--unique', required=True, help='File to find unique families from')
parser.add_argument('--file2', required=True, help='Second TSV file to compare against')
parser.add_argument('--file3', required=True, help='Third TSV file to compare against')
parser.add_argument('--output', required=True, help='Output file path')

args = parser.parse_args()
analyze_unique_families(args.unique, args.file2, args.file3, args.output)

args = parser.parse_args()
analyze_unique_families(args.unique, args.file2, args.file3, args.output)
