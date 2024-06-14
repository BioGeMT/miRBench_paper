# Precision-Recall Curve

This directory contains the code to generate the precision-recall curves for the miRNA target prediction methods.

## Requirements

- Python 3.x
- pandas
- matplotlib
- scikit-learn

## Usage

```bash
python pr_curves.py --ifile <input_file> --ofile <output_file> [--predictors <list_of_predictors> --title <plot_title>]
```

### Arguments

- `--ifile`: Input file containing the prediction scores in TSV format (default: STDIN)
- `--ofile`: Output file to save the precision-recall curves (default: STDOUT)
- `--predictors`: List of predictors to plot (default: all predictors)
- `--title`: Title for the plot (default: Precision-Recall Curve)

### Example

```bash
python pr_curves.py \
    --ifile pr_curves_example_input/Hejret_2023_miRNA_test_set_1_predictions.tsv \
    --ofile pr_curves_example_output/Hejret_2023_miRNA_test_set_1_pr_curves.png \
    --title "Hejret 2023 miRNA 1:1 test set"
```


