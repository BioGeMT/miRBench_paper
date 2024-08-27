# Negative Example Generator for miRNA-Gene Interactions

This Python script generates negative examples for miRNA-gene interaction datasets. 

## Features

- Generates negative examples based on positive miRNA-gene interaction data
- Ensures a minimum edit distance between positive and negative miRNA examples
- Supports flexible negative-to-positive example ratios
- Handles large datasets efficiently by processing data in blocks
- Preserves consistency in feature and test attributes within gene blocks

## Requirements

- Python 3.x
- pandas
- Levenshtein

You can install the required packages using pip:

`pip install pandas python-Levenshtein`

## Usage

Run the script from the command line with the following arguments:

python make_neg_sets.py --ifile INPUT_FILE --ofile OUTPUT_FILE [--neg_ratio RATIO] [--min_required_edit_distance DISTANCE]

Arguments:
- `--ifile`: Path to the input file (positive examples, must be sorted by 'gene')
- `--ofile`: Path to the output file
- `--neg_ratio`: Number of negative examples to generate per positive example (default: 1)
- `--min_required_edit_distance`: Minimum required edit distance for negative examples (default: 3)

## Input File Format

The input file should be a tab-separated file with the following columns:
- gene
- noncodingRNA
- noncodingRNA_fam
- feature
- test
- label

The file MUST be sorted by the 'gene' column.

## Output

The script generates an output file in the same format as the input, including both positive and negative examples. Negative examples are labeled with 0.

## Notes

- The script uses a fixed random seed (42) for reproducibility.
- Warnings are printed to stderr for any inconsistent blocks or if the requested number of negative examples couldn't be generated.
- Execution time is printed at the end of the process.
