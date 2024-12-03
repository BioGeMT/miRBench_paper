import pandas as pd
import sys
import argparse
import warnings

def filter_and_create_table(data):

    # filter rows where "noncodingRNA_type" is "miRNA".
    filtered_data = data[data['noncodingRNA_type'] == 'miRNA']

    # create the new dataframe with specific column names and transformations.
    filtered_table = pd.DataFrame({
        'gene': filtered_data['seq.g'],
        'noncodingRNA': filtered_data['noncodingRNA_seq'],
        'noncodingRNA_name': filtered_data['noncodingRNA'],
        'noncodingRNA_fam': filtered_data['noncodingRNA_fam'].apply(lambda x: x if x != '0' else 'unknown'),
        'feature': filtered_data['feature'],
        'test': filtered_data['chr.g'].apply(lambda x: True if x == '1' else False),
        'label': '1',
        'chr': filtered_data['chr.g'],
        'start': filtered_data['start.g'],
        'end': filtered_data['end.g'],
        'strand': filtered_data['strand.g']
    }, columns=['gene', 'noncodingRNA', 'noncodingRNA_name', 'noncodingRNA_fam', 'feature', 'test', 'label', 'chr', 'start', 'end', 'strand'])

    return filtered_table

# read input data from a file or stdin.

def read_input(input_file):

    try:
        if input_file:
            data = pd.read_csv(input_file, sep='\t')
        else:
            data = pd.read_csv(sys.stdin, sep='\t')
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
    
    parser = argparse.ArgumentParser(description="miRNA Finder")
    parser.add_argument('--ifile', help="Input file (default: STDIN)")
    parser.add_argument('--ofile', help="Output file (default: STDOUT)")

    args = parser.parse_args()

    # read input data.
    data = read_input(args.ifile)

    # process data.
    filtered_table = filter_and_create_table(data)

    # write output data.
    write_output(filtered_table, args.ofile)

if __name__ == "__main__":
    main()
