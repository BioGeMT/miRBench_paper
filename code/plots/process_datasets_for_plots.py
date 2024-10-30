import pandas as pd
import argparse
from datetime import datetime

def process_hejret(input_filepath):
    # Load the Hejret dataset
    hejret = pd.read_csv(input_filepath, sep='\t')
    # Filter for miRNA only
    hejret = hejret[hejret['noncodingRNA_type'] == 'miRNA']
    # Drop duplicates in gene + miRNA sequence but keep first occurrence
    hejret = hejret.drop_duplicates(subset=['seq.g', 'noncodingRNA_seq'], keep='first')
    hejret = hejret.rename(columns={'seq.g': 'gene',
                                    'noncodingRNA': 'noncodingRNA_renamed',
                                    'noncodingRNA_seq': 'noncodingRNA',
                                    'miRNA_fam': 'noncodingRNA_fam'})
    print("Hejret2023: ", len(hejret))
    return hejret

def process_klimentova_full(input_filepath):
    # Load the Klimentova dataset
    klimentova_full = pd.read_csv(input_filepath, sep='\t')
    klimentova_full = klimentova_full[klimentova_full['noncodingRNA_type'] == 'miRNA']
    klimentova_full = klimentova_full.drop_duplicates(subset=['seq.g', 'noncodingRNA_seq'], keep='first')
    klimentova_full = klimentova_full.rename(columns={'seq.g': 'gene',
                                                      'noncodingRNA': 'noncodingRNA_renamed',
                                                      'noncodingRNA_seq': 'noncodingRNA',
                                                      'miRNA_fam': 'noncodingRNA_fam'})
    return klimentova_full

def process_klimentova_test(input_filepath):
    # Load the Klimentova test dataset
    klimentova_test = pd.read_csv(input_filepath, sep='\t')
    klimentova_test = klimentova_test[klimentova_test['label'] == 1]
    return klimentova_test.rename(columns={'miRNA_fam': 'noncodingRNA_fam'})

def merge_datasets(klimentova_full, klimentova_test):
    # Merge Klimentova full and test datasets
    klimentova_final = klimentova_full.merge(klimentova_test, on=['gene', 'noncodingRNA'], how='inner')
    print("Klimentova2022: ", len(klimentova_final))
    return klimentova_final

def process_manakov(input_filepath):
    # Load the Manakov dataset
    manakov = pd.read_csv(input_filepath, sep='\t')
    # Filter for positive labels only
    manakov = manakov[manakov['label'] == 1]
    print("Manakov2022: ", len(manakov))
    return manakov

def main(hejret_input, klimentova_full_input, klimentova_test_input, manakov_input, hejret_output, klimentova_output, manakov_output):
    # Process Hejret dataset
    hejret = process_hejret(hejret_input)

    # Process and merge Klimentova datasets
    klimentova_full = process_klimentova_full(klimentova_full_input)
    klimentova_test = process_klimentova_test(klimentova_test_input)
    klimentova = merge_datasets(klimentova_full, klimentova_test)

    # Process Manakov dataset
    manakov = process_manakov(manakov_input)

    # Save to CSV files
    hejret.to_csv(hejret_output, sep='\t', index=False)
    klimentova.to_csv(klimentova_output, sep='\t', index=False)
    manakov.to_csv(manakov_output, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process miRNA datasets.')
    parser.add_argument('--hejret_input', type=str, required=True, help='Input file path for Hejret dataset')
    parser.add_argument('--klimentova_full_input', type=str, required=True, help='Input file path for Klimentova raw HD output dataset')
    parser.add_argument('--klimentova_test_input', type=str, required=True, help='Input file path for Klimentova model testing dataset')
    parser.add_argument('--manakov_input', type=str, required=True, help='Input file path for Manakov dataset')
    parser.add_argument('--hejret_output', type=str, required=True, help='Output file path for processed Hejret dataset')
    parser.add_argument('--klimentova_output', type=str, required=True, help='Output file path for processed Klimentova dataset')
    parser.add_argument('--manakov_output', type=str, required=True, help='Output file path for processed Manakov dataset')

    args = parser.parse_args()
    
    main(args.hejret_input, args.klimentova_full_input, args.klimentova_test_input, args.manakov_input,
         args.hejret_output, args.klimentova_output, args.manakov_output)

