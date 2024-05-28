
# Feature Prevalence Calculator

This script calculates the prevalence of features in a dataset, outputs the results to a TSV file, and optionally creates a bar plot of the feature prevalence.

## Requirements

- Python 3.x
- pandas
- matplotlib
- argparse

## Usage

```bash
python feature_prevalence.py --ifile <input_file> --pfile <plot_file> --ofile <output_tsv_file>
```

### Arguments

- `--ifile`: Input TSV file (default: stdin)
- `--pfile`: Output plot image file (default: none)
- `--ofile`: Output TSV file (default: stdout)

### Example

```bash
python feature_prevalence.py --ifile data.tsv --pfile feature_prevalence.png --ofile feature_prevalence.tsv
```

## Description

- The script reads the input data from a file or stdin.
- Calculates the prevalence of each feature in the dataset.
- Outputs the feature prevalence to a TSV file or stdout.
- Optionally creates a bar plot of the feature prevalence.

## Functions

- `calculate_feature_prevalence(data)`: Calculates the prevalence of each feature in the dataset.
- `read_input(input_file)`: Reads the input data from a file or stdin.
- `write_output(data, output_file)`: Writes the processed data to a TSV file or stdout.
- `plot_feature_prevalence(data, plot_file)`: Creates and saves a bar plot of the feature prevalence.
- `main()`: Main function to handle argument parsing and calling the other functions.
