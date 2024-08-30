## Description

The `venn.py` script processes TSV files containing miRNA family data (family_counts.tsv). It creates a Venn diagram to visualize the overlap of miRNA families between these studies, calculates various statistics, and outputs the results to files.

## Usage

Run the script from the command line with the following syntax:

```
python venn.py -i input_file.tsv -o output_diagram.png -s output_stats.tsv
```

### Arguments:
- `-i`, `--input`: Required. The input TSV file containing miRNA family data.
- `-o`, `--output`: Required. The name of the output file for the Venn diagram (e.g., diagram.png).
- `-s`, `--stats`: Required. The name of the output TSV file for the statistics.

## Input File Format

The input TSV file should have a header row and include a column named 'miRNA Family' which contains the miRNA family names. It should also have columns named 'Manakov2022', 'Hejret2023', and 'Klimentova2022' with their respective count values.

## Dependencies

- pandas
- matplotlib
- matplotlib_venn
- argparse

