
# Scatter Plot Creator

This script creates a scatter plot to visualize the relationship between the number of occurrences of each miRNA and the percentage of perfect seed interaction.

## Requirements

- Python 3.x
- pandas
- matplotlib
- numpy
- argparse

## Usage

```bash
python scatter_plot.py --ifile <input_file> --ofile <output_file>
```

### Arguments

- `--ifile`: Input TSV file (default: STDIN)
- `--ofile`: Output image file (default: STDOUT)

### Example

```bash
python scatter_plot.py --ifile data.tsv --ofile scatterplot.png
```

## Description

- The script reads the input data from a file or stdin.
- Transforms the number of occurrences using a logarithmic base 10 scale.
- Creates a scatter plot with the transformed number of occurrences on the x-axis and the percentage of perfect seed interaction on the y-axis.
- Saves the plot to a file or displays it.

## Functions

- `plot_data(input_file, output_file)`: Reads the data, processes it, and creates the scatter plot.
- `main()`: Main function to handle argument parsing and calling the `plot_data` function.
