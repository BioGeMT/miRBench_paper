# Evaluation of predictions from various predictors

## Requirements 

- scikit-learn
- pandas
- numpy

## Usage

```bash
python eval_metrics.py --ifile <input_predictions_file> --ofile <output_metrics_file> [--predictors <list_of_predictors>]`
```

### Arguments

- `--ifile`: Input file containing the prediction scores in TSV format (default: STDIN)
- `--ofile`: Output file to save the evaluation metrics (default: STDOUT)
- `--predictors`: List of predictors to evaluate (default: all predictors)

### Example

```bash
python eval_metrics.py \
--ifile example_input/AGO2_CLASH_Hejret2023_1_predictions.tsv \
--ofile example_output/AGO2_CLASH_Hejret2023_1_metrics.tsv
```
