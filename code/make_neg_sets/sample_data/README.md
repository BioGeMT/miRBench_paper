# Dummy Data Generator for miRNA-Gene Interactions

This Python script generates dummy data for miRNA-gene interaction datasets. It's designed to create synthetic datasets that mimic the structure of real miRNA-gene interaction data, which can be useful for testing, development, or demonstration purposes.

## Features

- Generates random RNA sequences for genes and miRNAs
- Creates a set of dummy miRNA families
- Assigns random features and test conditions
- Produces a specified number of positive examples
- Sorts the df based on the 'gene' column
- Outputs data in a tabular format suitable for generating negative examples

## Requirements

- Python 3.x
- pandas

## Usage

1. Set the desired number of samples in the `num_samples` variable.
2. Run the script: `python make_dummy_data.py`

The script will generate the dummy data and save it to a file named 'dummy_data_in.tsv' in the same directory.

## Output Format

The script generates a tab-separated file with the following columns:
- gene: Random DNA sequence of length 50
- noncodingRNA: Random DNA sequence of length 22 (miRNA)
- noncodingRNA_fam: miRNA family (e.g., 'mir-1', 'mir-2', etc.)
- feature: Random choice from ['three_prime_utr', 'five_prime_utr', 'intron', 'exon']
- test: Random boolean value (True or False)
- label: Always 1 (representing positive examples)

Example of output file is `dummy_data_in.tsv`.

## Notes

- This script generates only dummy positive examples (label = 1). 
- The output file is sorted by the 'gene' column in preparation for generating negatives with the `make_neg_sets.py` script. 