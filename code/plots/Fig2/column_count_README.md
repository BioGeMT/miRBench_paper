# Column Counter Script

This script processes multiple TSV files to count occurrences of values in a specified column. It analyzes miRNA data, supporting both family-based and sequence-based counting. Its output is necessary for venn.py and top_families.py scripts.

## Description

The script processes multiple TSV files containing miRNA data and counts the occurrences of values in a specified column (either miRNA families or sequences). It generates a summary table in TSV format, where rows represent unique values from the specified column and columns show their counts in each input file.

## Usage

Run the script from the command line with the following syntax:

```bash
python column_count.py input_file1.tsv input_file2.tsv input_file3.tsv output_file.tsv count_column
```

### Arguments
* `input_file1.tsv`, `input_file2.tsv`, `input_file3.tsv`: Three input TSV files containing miRNA data
* `output_file.tsv`: The name of the output TSV file where the results will be written
* `count_column`: Name of the column to count values from (e.g., 'noncodingRNA_fam' for families or 'noncodingRNA' for sequences)

### Common Use Cases

For counting miRNA families:
```bash
python column_count.py input1.tsv input2.tsv input3.tsv family_counts.tsv noncodingRNA_fam
```

For counting miRNA sequences:
```bash
python column_count.py input1.tsv input2.tsv input3.tsv sequence_counts.tsv noncodingRNA
```

## Input File Format

The input TSV files should:
* Have a header row
* Include the specified column (either 'noncodingRNA_fam' for families or 'noncodingRNA' for sequences)
* Be tab-separated

## Output

The script generates a TSV file with the following structure:
* First column contains the unique values from the specified count column
* Subsequent columns contain the count of each value in each input file
* Column headers are named according to input datasets: [count_column, input1, input2, input3]
* Rows are sorted by total frequency across all files (descending order)

## Dependencies
* Python 3.x
* Standard Python libraries: argparse, csv, collections, os
