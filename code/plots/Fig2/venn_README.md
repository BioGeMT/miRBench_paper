## Description

The `venn.py` script processes TSV files containing unique miRNA sequences data (sequence_counts.tsv). It creates a Venn diagram to visualize the overlap of unique miRNA sequences between these studies, calculates various statistics, and outputs the results to files.

## Usage

Run the script from the command line with the following syntax:

```
python venn.py -i input_file.tsv -o output_diagram.png -s output_stats.tsv -d detailed_stats.tsv
```

### Arguments:
- `-i`, `--input`: Required. The input TSV file containing unique miRNA sequences data.
- `-o`, `--output`: Required. The name of the output file path for the Venn diagram (e.g., diagram.png).
- `-s`, `--stats`: Required. The name of the output TSV file  path for the statistics.
- `-d`, `--detailed`: Required. The name of the output TSV file path with detailed intersection statistics. 

## Input File Format

The input TSV file should have a header row and include a column named 'noncodingRNA_sequence' which contains the unique miRNA sequences. It should also have columns named 'Manakov2022', 'Hejret2023', and 'Klimentova2022' with their respective count values.

## Dependencies

- pandas
- matplotlib
- matplotlib_venn
- argparse

