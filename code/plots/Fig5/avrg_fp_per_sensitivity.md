# Average Scanning Length per Sensitivity Threshold Plot

## Description
The `avrg_fp_per_sensitivity_threshold.py` script generates a plot that shows the average scanning length (1/FPR) for various prediction models at a specified sensitivity threshold.

## Usage
Run the script from the command line with the following syntax:
python avrg_fp_per_sensitivity_threshold.py --dataset-dir <dataset_directory> --dataset-name <dataset_name>

### Arguments:
- `--dataset-dir`: (Optional) The directory containing the dataset file. Default is `'predictions/'`.
- `--output-dir`: (Optional) The directory where the output figures and intermediate results will be saved. Default is `'output/'`.
- `--dataset-name`: (Optional) The name of the dataset file. Default is `'AGO2_eCLIP_Manakov2022_100_CNN_predictions'`.
- `--thresholds`: (Optional) A list of desired sensitivity thresholds (recall) to use for the plot. Default is `[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]`.
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
- json
- argparse
- random
- pandas
- numpy
- matplotlib
- scikit-learn (for `roc_curve`)

Make sure you have these libraries installed before running the script.

## Output
The script will generate two output files in the `'output/'` directory:
- `<dataset_name>.avrg_fp_per_sensitivity_threshold.svg`: A Scalable Vector Graphics (SVG) file containing the average scanning length per sensitivity threshold plot.
- `<dataset_name>.avrg_fp_per_sensitivity_threshold.png`: A Portable Network Graphics (PNG) file containing the average scanning length per sensitivity threshold plot.

Additionally, the script will save the scanning length results in a JSON file named `'scanning_length_results.json'` in the `'output/'` directory.