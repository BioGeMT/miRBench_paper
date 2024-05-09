import pandas as pd
import sys
import os

def filter_and_create_table(input_file, output_file='filtered_data.tsv'):
    
    # ensure the output file has a .tsv extension
    if not output_file.endswith('.tsv'):
        output_file = os.path.splitext(output_file)[0] + '.tsv'

    # load the input .tsv file
    data = pd.read_csv(input_file, sep='\t')
    
    # filter rows where "noncodingRNA_type" is "miRNA"
    filtered_data = data[data['noncodingRNA_type'] == 'miRNA']
    
    # create the new dataframe with some new column names
    filtered_table = pd.DataFrame({
        'id': range(1, len(filtered_data) + 1),
        'seq.g': filtered_data['seq.g'],
        'seq.m': filtered_data['noncodingRNA_real_seq'],
        'noncodingRNA_fam': filtered_data['noncodingRNA_fam'],
        'feature': filtered_data['feature'],
        'test': filtered_data['chr.g'].apply(lambda x: True if x == '1' else False),
        'label': '1'
    })
    
    # save the new dataframe to a .tsv file
    filtered_table.to_csv(output_file, sep='\t', index=False)
    print(f"Filtered table saved to {output_file}")

if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print("no input file")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else 'filtered_data.tsv'
        filter_and_create_table(input_file, output_file)
