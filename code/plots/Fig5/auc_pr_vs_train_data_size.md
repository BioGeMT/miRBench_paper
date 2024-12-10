# AUC-PR vs Dataset Size Plot

## Description
The `auc_pr_vs_dataset_size.py` script generates a plot that visualizes the relationship between the Area Under the Precision-Recall Curve (AUC-PR) and the dataset size for various prediction methods.

## Usage
Run the script from the command line with the following syntax:
python auc_pr_vs_dataset_size.py --dataset-dir <dataset_directory> --dataset-name <dataset_name>

### Arguments:
- `--dataset-dir`: (Optional) The directory containing the dataset file. Default is `'predictions/'`.
- `--dataset-name`: (Optional) The name of the dataset file. Default is `'AGO2_eCLIP_Manakov2022_100_CNN_predictions'`.
- `--manakov-rename-map`: (Optional) A dictionary that maps the original Manakov dataset column names to the desired names for the plot. Default is: {
'CNN_Manakov_full': 'Manakov_2,524,246',
"CNN_Manakov_subset_200": 'Manakov_subset_200',
"CNN_Manakov_subset_2k": 'Manakov_subset_2k',
"CNN_Manakov_subset_7720": 'Manakov_subset_7720',
"CNN_Manakov_subset_20k": 'Manakov_subset_20k',
"CNN_Manakov_subset_200k": 'Manakov_subset_200k',
}
- `--method-names`: (Optional) A list of method names to include in the plot. Default is: ["random", "Manakov_subset_200", "Manakov_subset_2k", "Manakov_subset_7720",
"Manakov_subset_20k", "Manakov_subset_200k", "Manakov_2,524,246"]

## Input File Format
The script expects a TSV file (tab-separated values) with a header row. The file should contain the following columns:
- `label`: The ground truth label column
- Additional columns for each prediction method, where the column name matches the values in the `--method-names` argument.

## Dependencies
The script requires the following Python libraries:
- random
- numpy
- pandas
- matplotlib
- scikit-learn (for `precision_recall_curve` and `auc`)

Make sure you have these libraries installed before running the script.

## Output
The script will generate two output files in the `'output/'` directory:
- `<dataset_name>.auc_pr_vs_train_data_size.svg`: A Scalable Vector Graphics (SVG) file containing the AUC-PR vs dataset size plot.
- `<dataset_name>.auc_pr_vs_train_data_size.png`: A Portable Network Graphics (PNG) file containing the AUC-PR vs dataset size plot.