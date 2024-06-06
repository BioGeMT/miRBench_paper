
# Box Plot Creator

This script creates a box plot to visualize the statistics of different noncoding RNA families and their percentage of perfect seed interaction, and saves the plot data to a TSV file.

## Requirements

- Python 3.x
- pandas
- matplotlib
- argparse

## Usage

```bash
python box_plot.py --ifile <input_file> --pfile <plot_file> --ofile <output_tsv_file> --min_interactions <min_interactions>
```

### Arguments

- `--ifile`: Input TSV file (default: stdin)
- `--pfile`: Output plot image file (default: save plot)
- `--ofile`: Output TSV file (default: stdout)
- `--min_interactions`: Minimum number of interactions to filter families (default: 10)

### Example

```bash
python box_plot.py --ifile data.tsv --pfile boxplot.png --ofile family_stats.tsv --min_interactions 10
```

## Description

- The script reads the input data from a file or stdin.
- Filters the families with a specified minimum number of interactions.
- Creates a box plot showing the mean, minimum, and maximum percentage of perfect seed interaction for each noncoding RNA family.
- Saves the plot to a file if specified.
- Saves the filtered family statistics to a TSV file, sorted by mean percentage in descending order.

## Functions

- `calculate_family_statistics(data, min_interactions)`: Filters the data by the minimum number of interactions and sorts it for plotting and TSV output.
- `read_input(input_file)`: Reads the input data from a file or stdin.
- `write_output(data, output_file)`: Writes the processed data to a TSV file or stdout.
- `plot_family_statistics(data, plot_file)`: Creates and saves the box plot.
- `main()`: Main function to handle argument parsing and calling the other functions.
