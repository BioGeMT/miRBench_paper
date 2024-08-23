# Evaluation of predictions from various tools

This directory contains code to generate evaluation metrics of predictions from various tools. 

## Requirements 

- scikit-learn
- pandas
- numpy

## For auc-pr or auc-roc

### Usage

```bash
python get_auc.py --ifile <input_predictions_file> --ofile <output_auc_file> [--predictors <list_of_predictors>] [--metric <auc_to_compute>]`
```

### Arguments

- `--ifile`: Input file containing the prediction scores in TSV format (default: STDIN)
- `--ofile`: Output file to save the aucs (default: STDOUT)
- `--predictors`: List of predictors to evaluate (default: all predictors)
- `--metric`: Evaluation metric to compute; auc-pr or auc-roc (default: auc-pr)

### Example

```bash
python get_auc.py \
--ifile example_input/AGO2_CLASH_Hejret2023_1_predictions.tsv \
--ofile example_output/AGO2_CLASH_Hejret2023_1_auc-pr.tsv
```

## For confusion matrix (cm) - TN, FP, FN, TP

### Usage

```bash
python confusion_matrix.py --ifile <input_predictions_file> --ofile <output_cm_file> [--predictors <list_of_predictors>]
```

### Arguments

- `--ifile`: Input file containing the prediction scores in TSV format (default: STDIN)
- `--ofile`: Output file to save the cms (default: STDOUT)
- `--predictors`: List of predictors to evaluate (default: all predictors)

### Example

```bash
python confusion_matrix.py \
--ifile example_input/AGO2_CLASH_Hejret2023_1_predictions.tsv \
--ofile example_output/AGO2_CLASH_Hejret2023_1_cm.tsv
```
