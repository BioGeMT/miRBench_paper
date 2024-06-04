
# Scatter Plot Creator

This script creates a scatter plot to visualize the relationship between the number of occurrences of each miRNA and the percentage of perfect seed interaction, and saves the plot data to a TSV file.

## Requirements

- Python 3.x
- pandas
- matplotlib
- numpy
- argparse

## Usage

```bash
python scatter_plot.py --ifile <input_file> --pfile <plot_file> --ofile <output_tsv_file>
```

### Arguments

- `--ifile`: Input TSV file (default: stdin)
- `--pfile`: Output plot image file (default: none)
- `--ofile`: Output TSV file (default: stdout)

### Example

```bash
python scatter_plot.py --ifile data.tsv --pfile scatterplot.png --ofile plot_data.tsv
```

## Description

- The script reads the input data from a file or stdin.
- Transforms the number of occurrences using a logarithmic base 10 scale.
- Creates a scatter plot with the transformed number of occurrences on the x-axis and the percentage of perfect seed interaction on the y-axis.
- Saves the plot to a file if specified.
- Saves the plot data to a TSV file, sorted by the number of occurrences in descending order.

## Functions

- `plot_data(input_file, plot_file, tsv_file)`: Reads the data, processes it, creates the scatter plot, and saves the plot data to a TSV file.
- `main()`: Main function to handle argument parsing and calling the `plot_data` function.
