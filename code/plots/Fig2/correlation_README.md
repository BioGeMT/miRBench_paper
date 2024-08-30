## Description

The `correlation.py` script processes multiple TSV files containing miRNA family data, calculates correlations between datasets, and generates a correlation heatmap. It also outputs the correlation data to a TSV file.

## Usage

Run the script from the command line with the following syntax:

```
python correlation.py input1.tsv input2.tsv input3.tsv output_correlations.tsv output_heatmap.png
```

### Arguments:
- `input1.tsv`, `input2.tsv`, `input3.tsv`: Three input TSV files containing miRNA family data.
- `output_correlations.tsv`: The name of the output TSV file where the correlation matrix will be written.
- `output_heatmap.png`: The name of the output PNG file for the correlation heatmap.

## Input File Format

The input TSV files should have a header row and include a column named 'miRNA_fam' which contains the miRNA family names.

## Output

The script generates two outputs:

1. A TSV file containing the correlation matrix of miRNA family percentages across the input datasets.
2. A PNG file with a heatmap visualization of the correlation matrix.

## Dependencies

The script requires the following Python libraries:
- pandas
- numpy
- seaborn
- matplotlib

