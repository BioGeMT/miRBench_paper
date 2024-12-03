## Description

The `seq_count.py` script processes multiple TSV files containing miRNA sequences data. It counts the occurrences of unique miRNA sequences in each file and generates a summary table. The output is a new TSV file with miRNA sequences as rows and input files as columns, showing the count of each miRNA sequence in each file. This output is necessary for venn.py

## Usage

Run the script from the command line with the following syntax:

```
python seq_count.py input_file1.tsv input_file2.tsv input_file3.tsv output_file.tsv
```

### Arguments:
- `input_file1.tsv`, `input_file2.tsv`, `input_file3.tsv`: Three input TSV files containing miRNA data.
- `output_file.tsv`: The name of the output TSV file where the results will be written.

## Input File Format

The input TSV files should have a header row and include a column named 'noncodingRNA' which contains the miRNA sequences.

## Output

The script generates a TSV file with the following structure:
- The first column contains unique miRNA sequences.
- Subsequent columns contain the count of each unique miRNA sequence in each input file.
- Column headers are derived from the input file names (without the .tsv extension).
- Rows are sorted based on the frequency of unique miRNA sequences in the first input file.
