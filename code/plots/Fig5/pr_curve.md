# Precision-Recall Plot

## Description
The `pr_curve.py` script generates a precision-recall (PR) curve plot for various prediction models using the provided dataset.

## Usage
Run the script from the command line with the following syntax:
python pr_curve.py --dataset-dir <dataset_directory> --dataset-name <dataset_name>
### Arguments:
- `--dataset-dir`: (Optional) The directory containing the dataset file. Default is `'predictions/'`.
- `--output-dir`: (Optional) The directory where the output figures and intermediate results will be saved. Default is `'output/'`.
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
- argparse
- json
- numpy
- pandas
- matplotlib
- scikit-learn (for `precision_recall_curve` and `precision_recall_fscore_support`)

Make sure you have these libraries installed before running the script.

## Output
The script will generate the following output files in the `'output/'` directory:
- `<dataset_name>.pr_curve.svg`: A Scalable Vector Graphics (SVG) file containing the precision-recall curve plot.
- `pr_curve_results.json`: A JSON file containing the area under the precision-recall curve (AUC-PR) values for each prediction method.