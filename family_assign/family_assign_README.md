
# miRNA Family Assignment Script

This script processes miRNA data to assign families to miRNA sequences based on mature sequence information. It reads input from specified files and outputs the results to a specified file.

## Requirements

- Python 3.x
- pandas

## Usage

```bash
python families.py --ifile <input_file> --mature <mature_file> --ofile <output_file>
```

### Arguments

- `--ifile`: Input file containing miRNA data in TSV format (default: STDIN)
- `--mature`: File containing mature miRNA sequences (default: STDIN) *downloaded from mirbase*
- `--ofile`: Output file to save the processed data (default: STDOUT)

### Example

```bash
python families.py --ifile data/filtered_data.tsv --mature data/mature.fa --ofile results/processed_output.tsv
```

## Description

- The script reads miRNA data from the input file and mature miRNA sequences from the mature file.
- It updates the `noncodingRNA_fam` column in the input data based on the mature sequences.
- It removes the 'hsa-' prefix from the `noncodingRNA_fam` values if present.
- The processed data is saved to the output file.

## Functions

- `read_input(file_path)`: Reads input data from a file using tab as a separator.
- `filter_and_create_table(data, mature_sequences)`: Processes the data to update the `noncodingRNA_fam` column based on mature sequences.
- `write_output(data, file_path)`: Writes the processed data to an output file using tab as a separator.
- `load_mature_sequences(file_path)`: Loads mature sequences from a file and maps sequences to their families.
- `main()`: Main function to handle argument parsing and calling the processing functions.
