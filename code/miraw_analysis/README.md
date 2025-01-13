# Pipeline for training models on k-mer count features from miRNA sequences

This pipeline consists of three scripts designed to be run after eachother, for training, prediction, and evaluation of models using k-mer count features extracted from miRNA sequences. 

The pipeline is designed to serve as an analysis of datasets to show there is a bias in miRNA distribution across classes in [miRAW datasets](https://bitbucket.org/bipous/miraw_data/src/master/) and in our group's datasets that predate the miRBench project. The bias is absent in the miRBench datasets as it has been mitigated against. 

miRAW datasets of interest include:

- miraw_data/PLOSComb/Data/DIANA_DataSet/DIANADataSet/
- miraw_data/PLOSComb/Data/ValidTargetSites/
- miraw_data/PLOSComb/Data/TestData/balanced10/

## Requirements

- Python 3.x
- Required libraries:
  - `pandas`
  - `numpy`
  - `scikit-learn`
  - `matplotlib`
  - `argparse`
  - `logging`
  - `pickle` (included in Python standard library)

Install the dependencies using:

```bash
pip install pandas numpy scikit-learn matplotlib
```

## Pipeline Overview

### 1. Training (train.py)
Trains a Decision Tree classifier using k-mer count features and saves the trained model.

#### Usage:
```bash
python train.py --train_set <path_to_training_dataset> --k <k_mer_length> --output_dir <output_directory>
```

#### Arguments:
- `--train_set`: Path to the training dataset (.tsv) containing miRNA sequences (column: `Mature_mirna_transcript`, `mature_miRNA_Transcript`, or `noncodingRNA`) and labels (column: `Positive_Negative`, `validation`, or `label`).
- `--k`: Length of k-mers to use for feature extraction.
- `--output_dir`: Directory to save the trained model and log file.

#### Output:
- Trained model saved as a `.pkl` file.
- Log file recording the training process.

---

### 2. Prediction (predict.py)
Loads trained models and generates predictions for a test dataset.

#### Usage:
```bash
python predict.py --test_set <path_to_test_dataset> --k <k_mer_length> --models <list_of_model_paths> --output_dir <output_directory>
```

#### Arguments:
- `--test_set`: Path to the test dataset (.tsv) containing miRNA sequences (column: `Mature_mirna_transcript`, `mature_miRNA_Transcript`, or `noncodingRNA`) and labels (column: `Positive_Negative`, `validation`, or `label`).
- `--k`: Length of k-mers to use for feature extraction.
- `--models`: List of paths to the trained model files (.pkl).
- `--output_dir`: Directory to save the predictions and log file.

#### Output:
- Test set `.tsv` file with additional columns containing predictions for each model and a random baseline.
- Log file recording the prediction process.

---

### 3. Evaluation (evaluate.py)
Evaluates the performance of predictions using Average Precision Score and Precision-Recall (PR) metrics, and generates PR curves.

#### Usage:
```bash
python evaluate.py --predictions <path_to_predictions_file> --output_dir <output_directory>
```

#### Arguments:
- `--predictions`: Path to output prediction `.tsv` files from the `predict.py` script, containing predictions for different models (column: `*_model`) and true labels (column: `Positive_Negative`, `validation`, or `label`).
- `--output_dir`: Directory to save PR curves and evaluation metrics.

#### Output:
- PR curve plots for each model saved as `.png` files.
- Evaluation metrics on each model (average precision score and PR AUC) saved as `.tsv` files.
- Log file recording the evaluation process.

## Example Pipeline Execution

1. Train a model:
```bash
python train.py --train_set ./example_input/AGO2_CLASH_Hejret2023_train.tsv --k 3 --output_dir ./example_output/models/
```

2. Predict using the trained model:
```bash
python predict.py --test_set ./example_input/AGO2_CLASH_Hejret2023_test.tsv --k 3 --models ./example_output/models/DecisionTree_3mer_AGO2_CLASH_Hejret2023_train.pkl --output_dir ./example_output/predictions/
```

3. Evaluate the predictions:
```bash
python evaluate.py --predictions ./example_output/predictions/AGO2_CLASH_Hejret2023_test_predictions.tsv --output_dir ./example_output/evaluation/
```
