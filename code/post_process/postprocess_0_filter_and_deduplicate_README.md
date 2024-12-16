# postprocess_0_filter_and_deduplicate_README

## Usage

Can be submitted to SLURM workload manager with `sbatch` or:
```bash
postprocess_0_filter_and_deduplicate.sh -i input_dir -o output_dir -n intermediate_dir
```
### Parameters

`-i`: Input directory containing .tsv files  
`-o`: Output directory for final files  
`-n`: Intermediate directory for intermediate files  

## Pipeline Steps

1. **Filtering**: Refer to ../filtering/README.md
2. **Deduplication**: Removes duplicates based on first two columns (gene and noncodingRNA sequence)

## Intermediate Files

- Filtered files: `{input_file}.filtered.tsv`

## Output Files

- Final files: `{input_file}.filtered.deduplicated.tsv`
- Log file: `postprocess_0_filter_and_deduplicate.log`