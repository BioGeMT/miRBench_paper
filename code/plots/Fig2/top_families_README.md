## Description

The `top_families.py` script processes TSV files containing miRNA family count data (family.counts.tsv). It calculates the percentage distribution of miRNA families, generates a bar plot to visualize the top 10 families across up to three datasets, and outputs the results to a TSV file.

## Usage

Run the script from the command line with the following syntax:

```
python top_families.py -i input_file.tsv -o output_plot.png -t output_data.tsv
```

### Arguments:
- `-i`, `--input`: Required. The input TSV file containing miRNA family counts.
- `-o`, `--output`: Required. The name of the output file for the plot (e.g., plot.png).
- `-t`, `--tsv`: Required. The name of the output TSV file for the plot results.

## Input File Format

The input TSV file should have a header row and include a column named 'miRNA Family' which contains the miRNA family names. Other columns should represent different datasets or conditions with their respective count values.

## Output

The script generates the following outputs:

1. A bar plot visualizing the distribution of the top 10 miRNA families across the input datasets.
   - The plot is saved to the specified output file.
   - Families are represented by grouped bars, color-coded for each dataset.

2. A TSV file containing:
   - Rows representing the top 10 miRNA families.
   - Columns for each dataset, showing both the raw count and calculated percentage for each family.

## Dependencies

- pandas
- matplotlib
- numpy
- argparse

