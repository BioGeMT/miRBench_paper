
# miRNA Seed Analysis

This script analyzes miRNA seed interactions and calculates summary statistics from TSV files. It can read input from a specified input folder or from standard input (stdin), and it outputs results to a specified output folder or to standard output (stdout).

## Requirements

- Python 3.x
- pandas
- biopython


## Usage

```bash
python seed_analysis.py --ifolder <input_folder> --ofolder <output_folder>
```

### Arguments

- `--ifolder`: Input folder containing TSV files (required if not using stdin)
- `--ofolder`: Output folder for results (required if not using stdout)

### Example

```bash
python seed_analysis.py --ifolder data/ --ofolder results/
```

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
