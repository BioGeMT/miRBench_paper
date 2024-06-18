# Pre-processing Pipeline

This repository contains a Bash script (`pre_process.sh`) that automates a data pre-processing pipeline. The pipeline includes several steps: filtering, family assignment, data splitting, and generating negative samples.

## Prerequisites

Before running the script, ensure you have the following installed:

- Bash
- Python 3
- Required Python scripts: `filtering.py`, `family_assign.py`, and `make_neg_sets.py`
- Necessary input files: `data.tsv` and `mature.fa`

## Directory Structure

The script assumes the following directory structure:

```
├── input_data
│   ├── data.tsv
│   ├── mature.fa
├── intermediate_data
├── output_data
├── filtering.py
├── family_assign.py
├── make_neg_sets.py
├── pre_process.sh
```

## Script Overview

The script performs the following steps:

1. **Setup Directories and Logging**: Creates `input_data`, `intermediate_data`, and `output_data` directories if they don't exist and sets up logging.

2. **Filtering**: Executes the filtering step using `filtering.py`.

3. **Family Assignment**: Runs the family assignment step using `family_assign.py`.

4. **Data Splitting**: Splits the data into training and testing sets based on a specified column.

5. **Negative Sample Generation**: Generates negative samples for both training and testing sets with different ratios using `make_neg_sets.py`.

## Usage

To run the pipeline, execute the following command:

```bash
bash pre_process.sh
```

### Step-by-step Process

1. **Filtering**:
    ```bash
    python3 filtering.py --ifile input_data/data.tsv --ofile intermediate_data/filtered_data.tsv
    ```

2. **Family Assignment**:
    ```bash
    python3 family_assign.py --ifile intermediate_data/filtered_data.tsv --mature input_data/mature.fa --ofile intermediate_data/family_assigned_data.tsv
    ```

3. **Data Splitting**:
    ```bash
    awk -F'	' 'NR==1{header=$0; print header > "intermediate_data/train_data.tsv"; print header > "intermediate_data/test_data.tsv"} NR>1{if($5=="False"){print > "intermediate_data/train_data.tsv"} else {print > "intermediate_data/test_data.tsv"}}' intermediate_data/family_assigned_data.tsv
    ```

4. **Negative Sample Generation**:
    ```bash
    for ratio in 1 10 100; do
        python3 make_neg_sets.py --ifile intermediate_data/train_data.tsv --ofile output_data/train_data_with_negatives_$ratio.tsv --neg_ratio $ratio
        python3 make_neg_sets.py --ifile intermediate_data/test_data.tsv --ofile output_data/test_data_with_negatives_$ratio.tsv --neg_ratio $ratio
    done
    ```

## Log File

The script logs all its output to `output_data/pipeline.log`.

## Conclusion

The `pre_process.sh` script automates a series of pre-processing tasks essential for preparing data for downstream analysis. Ensure all required files and scripts are in place before running the pipeline.
