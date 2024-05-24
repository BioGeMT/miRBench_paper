
# miRNA Finder

This script filters rows where `noncodingRNA_type` is `miRNA` from the input data and creates a new table with specified columns and transformations.

## Requirements

- Python 3.x
- pandas

## Installation

Install the required Python packages using pip:

```bash
pip install pandas
```

## Usage

```bash
python filtering.py --ifile <input_file> --ofile <output_file>
```

### Arguments

- `--ifile`: Input file (default: STDIN)
- `--ofile`: Output file (default: STDOUT)

### Example

```bash
python filtering.py --ifile data.tsv --ofile filtered_data.tsv
```

## Description

- The script reads the input data from a file or stdin.
- Filters rows where `noncodingRNA_type` is `miRNA`.
- Creates a new DataFrame with specific columns and transformations.
- Writes the output data to a file or stdout.

## Functions

- `read_input(input_file)`: Reads the input data.
- `filter_and_create_table(data)`: Filters and creates the new table.
- `write_output(data, output_file)`: Writes the output data.
- `main()`: Main function to handle argument parsing and calling the processing functions.
