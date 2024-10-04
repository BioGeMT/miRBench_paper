import pandas as pd
import numpy as np
import ast
import argparse

def is_valid_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def check_conservation_scores(x):
    try:
        scores = ast.literal_eval(x)
        return all(is_valid_number(score) for score in scores)
    except:
        return False

def check_conservation_length(row, column):
    try:
        conservation_scores = ast.literal_eval(row[column])
        target_length = len(row['gene'])
        return len(conservation_scores) == target_length and all(is_valid_number(score) for score in conservation_scores)
    except:
        return False

# clean data. 
def clean_file(df):
    # print(f"Initial number of rows: {len(df)}")

    # Drop rows with invalid 'target_phyloP' conservation scores
    df['valid_target_phyloP'] = df['target_phyloP'].apply(check_conservation_scores)
    invalid_phyloP_rows = df[~df['valid_target_phyloP']]
    df = df[df['valid_target_phyloP']]
    df = df.drop(columns=['valid_target_phyloP'])

    # Drop rows with invalid 'target_phastCons' conservation scores
    df['valid_target_phastCons'] = df['target_phastCons'].apply(check_conservation_scores)
    invalid_phastCons_rows = df[~df['valid_target_phastCons']]
    df = df[df['valid_target_phastCons']]
    df = df.drop(columns=['valid_target_phastCons'])

    print(f"Number of rows dropped due to invalid phyloP or phastCons conservation scores: {len(invalid_phyloP_rows) + len(invalid_phastCons_rows)}")

    # print(f"Number of rows remaining: {len(df)}")
    
    # # Check for NaN values in each column
    # print("\nNaN values per column:")
    # for column in df.columns:
    #     nan_count = df[column].isna().sum()
    #     print(f"{column}: {nan_count}")

    # Check if 'target_phyloP' scores match target length
    df['phyloP_check'] = df.apply(lambda row: check_conservation_length(row, 'target_phyloP'), axis=1)
    mismatched_phyloP_rows = df[~df['phyloP_check']]
    print(f"\nNumber of rows where 'target_phyloP' scores don't match target length: {len(mismatched_phyloP_rows)}")

    # Check if 'target_phastCons' scores match target length
    df['phastCons_check'] = df.apply(lambda row: check_conservation_length(row, 'target_phastCons'), axis=1)
    mismatched_phastCons_rows = df[~df['phastCons_check']]
    print(f"\nNumber of rows where 'target_phastCons' scores don't match target length: {len(mismatched_phastCons_rows)}")

    # Print a sample of mismatched rows
    if not mismatched_phyloP_rows.empty:
        print("\nSample of rows where 'target_phyloP' scores don't match target length:")
        print(mismatched_phyloP_rows.head())
        
        # Print details of the first mismatched row
        first_mismatch = mismatched_phyloP_rows.iloc[0]
        print("\nDetails of first mismatched row (phyloP):")
        print(f"Target length: {len(first_mismatch['gene'])}")
        try:
            phyloP_scores = ast.literal_eval(first_mismatch['target_phyloP'])
            print(f"'target_phyloP' scores length: {len(phyloP_scores)}")
        except:
            print("Unable to parse 'target_phyloP' scores")
        print(f"Target sequence: {first_mismatch['gene']}")
        print(f"'target_phyloP' scores: {first_mismatch['target_phyloP']}")

    if not mismatched_phastCons_rows.empty:
        print("\nSample of rows where 'target_phastCons' scores don't match target length:")
        print(mismatched_phastCons_rows.head())
        
        # Print details of the first mismatched row
        first_mismatch = mismatched_phastCons_rows.iloc[0]
        print("\nDetails of first mismatched row (phastCons):")
        print(f"Target length: {len(first_mismatch['gene'])}")
        try:
            phastCons_scores = ast.literal_eval(first_mismatch['target_phastCons'])
            print(f"'target_phastCons' scores length: {len(phastCons_scores)}")
        except:
            print("Unable to parse 'target_phastCons' scores")
        print(f"Target sequence: {first_mismatch['gene']}")
        print(f"'target_phastCons' scores: {first_mismatch['target_phastCons']}")

    # Drop the 'phyloP_check' and 'phastCons_check' columns
    df = df.drop(columns=['phyloP_check', 'phastCons_check'])

    # # Print final number of entries
    # print(f"Final number of rows: {len(df)}")

    return df

# read input data from a file or stdin.
def read_input(input_file):

    try:
        if input_file:
            data = pd.read_csv(input_file, sep='\t', low_memory=False)
        else:
            data = pd.read_csv(sys.stdin, sep='\t', low_memory=False)
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

    parser = argparse.ArgumentParser(description="Clean file after conservation scores have been added.")
    parser.add_argument('--ifile', help="Input file (default: STDIN)")
    parser.add_argument('--ofile', help="Output file (default: STDOUT)")

    args = parser.parse_args()

    # read input data.
    data = read_input(args.ifile)

    # process data.
    cleaned_file = clean_file(data)

    # write output data.
    write_output(cleaned_file, args.ofile)

if __name__ == "__main__":
    main()