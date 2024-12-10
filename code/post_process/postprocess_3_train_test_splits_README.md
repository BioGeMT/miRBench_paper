# postprocess_3_train_test_splits_README

## Usage

Can be submitted to SLURM workload manager with `sbatch` or:
```bash
postprocess_3_train_test_splits.sh -i input_dir -o output_dir
```
### Parameters

`-i`: Input directory containing .tsv files
`-o`: Output directory for final files

## Pipeline Steps

Splits TSV files into train and test sets based on 'test' column (6th column). Processes all files with '.negatives.tsv' extension, preserves headers, and outputs separate train/test files with appropriate suffixes.

## Output Files

- Final files: 
    - `{input_file}.train.tsv`
    - `{input_file}.test.tsv`
- Log file: `postprocess_3_train_test_splits.log`