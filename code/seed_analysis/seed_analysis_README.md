# miRNA Seed Analysis

This script analyzes miRNA seed interactions and calculates summary statistics from TSV files. It reads input from a specified input folder and outputs results to a specified output folder.

## Requirements

- Python 3.x
- pandas
- biopython

```bash
python seed_analysis.py <input_folder> <output_folder>

## Description

- The script loads all TSV files from the input folder or reads from stdin.
- Analyzes miRNA seed interactions and calculates perfect matches.
- Aggregates results and calculates summary statistics.
- Saves the aggregated results and summary statistics to the output folder or stdout.

## Functions

- `load_tsv_files(input_stream)`: Loads all TSV files from the input stream.
- `reverse_complement(seq)`: Returns the reverse complement of a given sequence.
- `perfect_match(mirna, target)`: Checks if the seed region of the microRNA matches the reverse complement of the target sequence.
- `process_file(df)`: Processes each file to get interaction counts and perfect matches.
- `aggregate_results(results)`: Aggregates results and calculates summary statistics.
- `save_results(aggregate_df, summary_stats, output_stream)`: Saves the results to the output stream.
- `main()`: Main function to handle argument parsing and calling the processing functions.
