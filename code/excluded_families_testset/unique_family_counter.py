import pandas as pd
import argparse

def analyze_unique_families(unique_input_file, file2_path, file3_path, output_path):
   df1 = pd.read_csv(unique_input_file, sep='\t')
   df2 = pd.read_csv(file2_path, sep='\t')
   df3 = pd.read_csv(file3_path, sep='\t')
   
   df1 = df1[(df1['noncodingRNA_fam'] != 'unknown') & (df1['noncodingRNA_fam'] != '0')]
   df2 = df2[(df2['noncodingRNA_fam'] != 'unknown') & (df2['noncodingRNA_fam'] != '0')]
   df3 = df3[(df3['noncodingRNA_fam'] != 'unknown') & (df3['noncodingRNA_fam'] != '0')]
   
   families1 = set(df1['noncodingRNA_fam'].unique())
   families2 = set(df2['noncodingRNA_fam'].unique())
   families3 = set(df3['noncodingRNA_fam'].unique())
   
   unique_to_1 = families1 - (families2 | families3)
   unique_counts = df1[df1['noncodingRNA_fam'].isin(unique_to_1)]['noncodingRNA_fam'].value_counts()
   
   with open(output_path, 'w') as f:
       f.write("noncodingRNA_fam\tcount\n")
       for family, count in unique_counts.items():
           f.write(f"{family}\t{count}\n")

   print(f"Unique families: {len(unique_to_1)}")
   print(f"Total occurrences: {unique_counts.sum()}")

parser = argparse.ArgumentParser()
parser.add_argument('--unique', required=True, help='File to find unique families from')
parser.add_argument('--file2', required=True, help='Second TSV file to compare against')
parser.add_argument('--file3', required=True, help='Third TSV file to compare against')
parser.add_argument('--output', required=True, help='Output file path')

args = parser.parse_args()
analyze_unique_families(args.unique, args.file2, args.file3, args.output)