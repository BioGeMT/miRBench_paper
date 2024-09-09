
## Description

The `families.py` script analyzes miRNA family prevalence and seed type distribution across multiple datasets. It processes the input data, calculates family counts and percentages, and generates both a summary TSV file and a plot visualizing the seed type distributions for the top miRNA families.

## Usage

Run the script from the command line with the following syntax:

```
python families.py dataset1.tsv dataset2.tsv dataset3.tsv -o output_plot.png -t output_summary.tsv
```

### Arguments:
- `dataset1.tsv`, `dataset2.tsv`, `dataset3.tsv`: Input TSV files containing miRNA data.
- `-o`, `--output`: (Optional) The name of the output PNG file for the plot. Default is 'seed_prevalence_plot.png'.
- `-t`, `--tsv_output`: (Optional) The name of the output TSV file for the summary results. Default is 'miRNA_family_analysis.tsv'.

## Input File Format

The input TSV files should have a header row and include columns named 'seq.g' (guide sequence), 'seq.m' (miRNA sequence), and 'miRNA_fam' (miRNA family).

## Output

The script generates two output files:
1. A TSV file containing counts and percentages of miRNA families across all datasets.
2. A PNG file with a grouped bar plot showing the distribution of seed types for the top 10 miRNA families.

## Dependencies

- pandas
- matplotlib
- argparse
- Bio (Biopython)

Make sure to have `seed_utils.py` in the same directory as this script.
