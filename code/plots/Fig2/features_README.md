## Description

The `features.py` script processes TSV files containing miRNA feature data. It categorizes the features, calculates their distribution percentages, and generates a stacked bar plot to visualize the distribution across up to three datasets. The script also can output the results to a TSV file.

## Usage

Run the script from the command line with the following syntax:

```
python features.py [input_file1.tsv] [input_file2.tsv] [input_file3.tsv] [-o output_plot.png] [-t output_data.tsv]
```

### Arguments:
- `[input_file1.tsv]`, `[input_file2.tsv]`, `[input_file3.tsv]`: Up to three input TSV files containing miRNA feature data.
- `-o`, `--output`: Optional. The name of the output file for the plot (e.g., plot.png).
- `-t`, `--tsv`: Optional. The name of the output TSV file for the calculated percentages.

## Input File Format

The input TSV files should have a header row and include a column named 'feature' which contains the miRNA feature names.

## Output

The script generates the following outputs:

1. A stacked bar plot visualizing the distribution of miRNA features across the input datasets.
   - The plot is displayed if no output file is specified, or saved to the specified file.
   - Features are color-coded and include: 5' UTR, 3' UTR, Intron, Exon, and Unknown.

2. A TSV file (if specified with `-t` option) containing:
   - Rows representing each input dataset.
   - Columns for each feature category with their respective percentages.

## Dependencies

- pandas
- matplotlib
- argparse


