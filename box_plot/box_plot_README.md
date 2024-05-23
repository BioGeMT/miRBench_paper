
# Box Plot Creator

This script creates a custom box plot for noncoding RNA families with mean interactions of 10 or more from the input data.

## Requirements

- Python 3.x
- pandas
- matplotlib
- seaborn


## Usage

```bash
python box_plot.py --ifile <input_file> --ofile <output_file>
```

### Arguments

- `--ifile`: Input file (default: STDIN)
- `--ofile`: Output file (default: STDOUT)
- `--min_imteractions`: Minimum interactions filter 

### Example

```bash
python box_plot.py --ifile data.tsv --ofile boxplot.png
```

## Description

- The script reads the input data from a file or stdin.
- Filters families with 10 or more interactions.
- Creates a custom box plot and saves it to a file or displays it.

## Functions

- `load_data(input_file)`: Loads the data from the input file.
- `add_lines(df, y_min, y_max, y_mean, width, color_min, color_max, color_mean)`: Adds lines for min, max, and mean values.
- `create_box_plot(data, output_file)`: Creates the box plot.
- `main()`: Main function to handle argument parsing and calling the processing functions.
