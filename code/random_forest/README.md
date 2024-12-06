
# Random Forest Model with k-mers

This script trains and evaluates a Random Forest model for classification tasks using k-mers generated from miRNA sequences. It plots the Precision-Recall (PR) curve to assess model performance. It is made to be trained and tested on the miRAW dataset (DIANADataset)

## Features
- Loads miRNA sequences ONLY, and associated binary labels (0/1).
- Generates k-mers from sequences for feature extraction.
- Trains a Random Forest classifier.
- Evaluates the model using AUCPR (Area Under the Precision-Recall Curve).
- Saves the PR curve to a specified output file.

## Prerequisites
- Python 3.6+
- Required Python packages:
  - `pandas`
  - `scikit-learn`
  - `matplotlib`

Install dependencies:
```bash
pip install pandas scikit-learn matplotlib
```

## Usage
### Command-line Arguments
- `--train`: Path to the training dataset file (required).
- `--test`: Path to the testing dataset file (required).
- `--output`: Path to save the Precision-Recall curve image (required).

### Example
```bash
python script_name.py --train train_data.tsv --test test_data.tsv --output PR_curve.png
```

### Input Format - DIANADataset
- Tab-separated `.tsv` file with two columns:
  - `Mature_mirna_transcript`: Sequence data (strings).
  - `Positive_Negative`: Binary labels (0/1).

### Output
- Prints the AUCPR value to the console.
- Saves the PR curve to the specified output file.
