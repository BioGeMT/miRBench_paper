
## Description

The `violin.py` script processes miRNA data from multiple datasets and generates a violin plot to visualize the prevalence of canonical seed types of unique miRNA families across datasets. It uses the seed matching logic defined in `seed_utils.py` to categorize miRNA-target interactions. The script now also outputs a TSV file containing the data used to create the violin plot.

## Usage

Run the script from the command line with the following syntax:

```
python violin.py dataset1.tsv dataset2.tsv dataset3.tsv -o output_plot.png -t output_data.tsv
```

### Arguments:
- `dataset1.tsv`, `dataset2.tsv`, `dataset3.tsv`: Input TSV files containing miRNA data.
- `-o`, `--output`: (Optional) The name of the output PNG file for the violin plot. Default is 'seed_prevalence_violin_plot.png'.
- `-t`, `--tsv`: (Optional) The name of the output TSV file for the processed data. Default is 'seed_prevalence_data.tsv'.

## Input File Format

The input TSV files should have a header row and include columns named 'gene' (guide sequence) and 'noncodingRNA' (miRNA sequence).


## Dependencies

- pandas
- matplotlib
- numpy
- argparse
- collections

Make sure to have `seed_utils.py` in the same directory as this script.

