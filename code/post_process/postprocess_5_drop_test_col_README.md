# postprocess_4_drop_test_col_README

## Usage

Can be submitted to SLURM workload manager with `sbatch` or:
```bash
postprocess_4_drop_test_col.sh -i input_dir -o output_dir
```
### Parameters

`-i`: Input directory containing .tsv files  
`-o`: Output directory for final files  

## Pipeline Steps

Removes the 'test' column (6th column).

## Output Files

- Final file: `{input_file}.drop_test_col.tsv`
- Log file: `postprocess_4_drop_test_col.log`