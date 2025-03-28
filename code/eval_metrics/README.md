# Evaluation of predictions from various tools

This directory contains code to generate evaluation metrics of predictions from various tools. 

## Requirements 

- scikit-learn
- pandas
- numpy

## Usage

```bash
python get_metric.py --ifile <input_predictions_file> --ofile <output_auc_file> --metric <metric_to_compute> [--predictors <list_of_predictors>]`
```

### Arguments

- `--ifile`: Input file containing the prediction scores in TSV format (default: STDIN)
- `--ofile`: Output file to save the aucs (default: STDOUT)
- `--metric`: Evaluation metric to compute; auc_pr, auc_roc or avg_p_score
- `--predictors`: List of predictors to evaluate (default: all predictors)

### Example

```bash
python get_metric.py \
--ifile example_input/AGO2_eCLIP_Klimentova2022_test_predictions.tsv \
--ofile example_output/AGO2_eCLIP_Klimentova2022_test_auc_pr.tsv
--metric auc_pr
```

