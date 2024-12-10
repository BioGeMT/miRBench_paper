# postprocess_4_add_conservation_README

## Usage

Can be submitted to SLURM workload manager with `sbatch` or:
```bash
postprocess_4_add_conservation.sh -i input_dir -o output_dir -p phyloP_file -c phastCons_file
```
### Parameters

`-i`: Input directory containing .tsv files
`-o`: Output directory for final files
`-p`: Path to phyloP file
`-c`: Path to phastCons file

## Pipeline Steps

Refer to ../conservation/README.md.

## Output Files

- Final files: `{input_file}.conservation.tsv`
- Log file: `postprocess_4_add_conservation.log`