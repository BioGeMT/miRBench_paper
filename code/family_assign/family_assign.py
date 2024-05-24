import pandas as pd
import argparse

def read_input(file_path):
    # read input data from a file using tab as a separator
    return pd.read_csv(file_path, sep='\t')

def filter_and_create_table(data, mature_sequences):
    # ensure the necessary column is present in the data
    if 'noncodingRNA_fam' not in data.columns:
        raise KeyError("The column 'noncodingRNA_fam' is not found in the input data.")
    
    # update 'noncodingRNA_fam' based on mature sequences if it is '0'
    data['noncodingRNA_fam'] = data.apply(
        lambda row: mature_sequences.get(row['seq.m'].replace('T', 'U'), row['noncodingRNA_fam']) if row['noncodingRNA_fam'] == '0' else row['noncodingRNA_fam'],
        axis=1
    )
    
    # remove 'hsa-' prefix from 'noncodingRNA_fam' values if present
    data['noncodingRNA_fam'] = data['noncodingRNA_fam'].apply(lambda x: x.replace('hsa-', '') if 'hsa-' in x else x)
    
    return data

def write_output(data, file_path):
    # write the processed data to an output file using tab as a separator
    data.to_csv(file_path, sep='\t', index=False)

def load_mature_sequences(file_path):
    # load mature sequences from a file and map sequences to their families
    mature_sequences = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            header = lines[i].strip()
            sequence = lines[i+1].strip().replace('T', 'U')
            family = header.split()[0][1:]
            mature_sequences[sequence] = family
    return mature_sequences

def main():
    # parse command-line arguments for input, mature, and output files
    parser = argparse.ArgumentParser(description="miRNA Finder")
    parser.add_argument('--ifile', help="Input file (default: STDIN)")
    parser.add_argument('--mature', help="Mature file (default: STDIN)")
    parser.add_argument('--ofile', help="Output file (default: STDOUT)")

    args = parser.parse_args()

    # load mature sequences from the specified file
    mature_sequences = load_mature_sequences(args.mature)

    # read input data from the specified file
    data = read_input(args.ifile)

    # process the data to update the 'noncodingRNA_fam' column
    filtered_table = filter_and_create_table(data, mature_sequences)

    # write the processed data to the specified output file
    write_output(filtered_table, args.ofile)

if __name__ == "__main__":
    main()
