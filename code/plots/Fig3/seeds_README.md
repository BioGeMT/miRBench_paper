
## Description

The `seeds.py` script analyzes miRNA seed prevalence across multiple datasets. It processes the input data, categorizes miRNA-target interactions based on seed types, and generates both a summary TSV file and a stacked bar plot visualizing the seed type distributions.

## Usage

Run the script from the command line with the following syntax:

```
python seeds.py dataset1.tsv dataset2.tsv dataset3.tsv -o output_plot.png -t output_summary.tsv
```

### Arguments:
- `dataset1.tsv`, `dataset2.tsv`, `dataset3.tsv`: Input TSV files containing miRNA data.
- `-o`, `--output`: (Optional) The name of the output PNG file for the plot. Default is 'seed_prevalence_plot.png'.
- `-t`, `--tsv`: (Optional) The name of the output TSV file for the summary results. Default is 'seed_prevalence_results.tsv'.

## Input File Format

The input TSV files should have a header row and include columns named 'gene' (guide sequence) and 'noncodingRNA' (miRNA sequence).

## Output

The script generates two output files:
1. A TSV file containing counts and percentages of each seed type for each dataset.
2. A PNG file with a stacked bar plot showing the distribution of seed types across datasets.

## Dependencies

- pandas
- matplotlib
- argparse
- collections

Make sure to have `seed_utils.py` in the same directory as this script.
