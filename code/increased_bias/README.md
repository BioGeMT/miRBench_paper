# Pipeline for evaluating tools on dataset with increased miRNA frequency bias

## Requirements

- Python 3.x
- Required libraries:
  - `pandas`
  - `numpy`
  - `scikit-learn`
  - `miRBench`
  - `biopython`
  - `viennarna`
  - `torch`
  - `tensorflow`
  - `typing-extension`

Install the dependencies using:

```bash
pip miRBench
pip install numpy==1.24.3 biopython==1.83 viennarna==2.7.0 torch==1.9.0 tensorflow==2.13.1 typing-extensions==4.5.0
pip install scikit-learn
```

## 1. Generate dataset with increased bias

This tools take `train` and `test` datasets and modifies test dataset to increase miRNA frequency bias. Only those positive samples are kept where miRNA in the sample is more frequent in the positive class in training samples. The same for negative samples.

```bash
python make_biased_dataset.py --input_train_file example_input/miRNA_train_set.tsv --input_test_file example_input/miRNA_test_set_1.tsv --output_file example_output/biased_test.tsv
```

## 2. Get predictions from tools in miRBench package

This tool runs all predictors in the miRBench package to get predictions on a provided test file.

```bash
python predict.py --test example_output/biased_test.tsv --output example_output/biased_test_predictions.tsv
```

## 3. Compute average precision score

To get evaluation metrics, we can call `get_metric.py` tool from `../eval_metrics`.

```bash
python ../eval_metrics/get_metric.py --ifile example_output/biased_test_predictions.tsv --ofile example_output/biased_test_avg_p_score.tsv --metric avg_p_score
```