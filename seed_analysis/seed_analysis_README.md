
# miRNA Seed Analysis

This script analyzes miRNA seed interactions and calculates summary statistics from multiple TSV files in a specified input folder.

## Requirements

- Python 3.x
- pandas
- biopython


## Usage

```bash
python seed_analysis.py --ifolder <input_folder> --ofolder <output_folder>
```

### Arguments

- `--ifolder`: Input folder containing TSV files (required)
- `--ofolder`: Output folder for results (required)

### Example

```bash
python seed_analysis.py --ifolder data/ --ofolder results/
```

## Description

- The script loads all TSV files from the input folder.
- Analyzes miRNA seed interactions and calculates perfect matches.
- Aggregates results and calculates summary statistics.
- Saves the aggregated results and summary statistics to the output folder.

## Functions

- `load_tsv_files(input_folder)`: Loads all TSV files from the input folder.
- `reverse_complement(seq)`: Returns the reverse complement of a given sequence.
- `perfect_match(mirna, target)`: Checks if the seed region of the microRNA matches the reverse complement of the target sequence.
- `process_file(df)`: Processes each file to get interaction counts and perfect matches.
- `aggregate_results(results)`: Aggregates results and calculates summary statistics.
- `save_results(aggregate_df, summary_stats, output_folder)`: Saves the results to the output folder.
- `main()`: Main function to handle argument parsing and calling the processing functions.
