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

def validate_conservation_scores(df):
    # Validate conservation score values
    df['target_phyloP'] = df['target_phyloP'].apply(lambda x: x if check_conservation_scores(x) else np.nan)
    df['target_phastCons'] = df['target_phastCons'].apply(lambda x: x if check_conservation_scores(x) else np.nan)
    # Validate conservation scores length
    df['target_phyloP'] = df.apply(lambda row: row['target_phyloP'] if check_conservation_length(row, 'target_phyloP') else np.nan, axis=1)
    df['target_phastCons'] = df.apply(lambda row: row['target_phastCons'] if check_conservation_length(row, 'target_phastCons') else np.nan, axis=1)
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


def main():
    parser = argparse.ArgumentParser(description="Validate conservation score values and length.")
    parser.add_argument('--ifile', help="Input file (default: STDIN)")
    parser.add_argument('--ofile', help="Output file (default: STDOUT)")

    # Parse arguments
    args = parser.parse_args()

    # Read input data into a df
    data = read_input(args.ifile)

    # Validate conservation score values and length
    validated_data = validate_conservation_scores(data)

    # Write output data
    write_output(validated_data, args.ofile)

if __name__ == "__main__":
    main()