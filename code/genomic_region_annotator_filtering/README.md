# Genomic Region Annotator Filtering

This script filters the `*.annotated_site_summary.tsv` output from genomic region annotation and creates a new table with the selected downstream columns in a fixed order.

## Requirements

- Python 3.x
- pandas


## Usage

```bash
python genomic_region_annotator_filtering.py --ifile <input_file> --ofile <output_file>
```

### Arguments

- `--ifile`: Input file (default: STDIN)
- `--ofile`: Output file (default: STDOUT)

### Example

```bash
python genomic_region_annotator_filtering.py --ifile sample.annotated_site_summary.tsv --ofile sample.annotated.filtered.tsv
```

## Description

- The script reads the input data from a file or stdin.
- Selects the required columns by name from the annotated site summary table.
- Renames the final four selected-transcript region columns for downstream use.
- Preserves only the requested output columns in the defined order.
- Writes the output data to a file or stdout.

## Output Columns

The script keeps the following columns in this order:

- `gene`
- `noncodingRNA`
- `noncodingRNA_name`
- `noncodingRNA_fam`
- `feature`
- `test`
- `label`
- `chr`
- `start`
- `end`
- `strand`
- `Nunique`
- `dominant_region`
- `regions_present`
- `read_start_in_sel_tx_1based`
- `read_end_in_sel_tx_1based`

The final four columns are renamed from:

- `dominant_region_selected`
- `regions_present_selected`
- `selected_read_start_in_tx_1based`
- `selected_read_end_in_tx_1based`

## Functions

- `read_input(input_file)`: Reads the input data.
- `filter_columns(data)`: Selects the required columns and renames the final four.
- `write_output(data, output_file)`: Writes the output data.
- `main()`: Main function to handle argument parsing and calling the processing functions.
