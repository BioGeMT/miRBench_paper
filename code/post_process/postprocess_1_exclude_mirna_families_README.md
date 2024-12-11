# postprocess_1_exclude_mirna_families_README

## Usage

Can be submitted to SLURM workload manager with `sbatch` or:
```bash
postprocess_1_exclude_mirna_families.sh -i primary_input_file -j input_file_2 -k input_file_3 -o output_dir -n intermediate_dir
```
### Parameters

`-i`: Primary input .tsv file to be split into an excluded (miRNA fams unique to this dataset) and a remaining dataset  
`-j`: First .tsv file to compare miRNA fams against  
`-k`: Second .tsv file to compare miRNA fams against  
`-o`: Output directory for final files  
`-n`: Intermediate directory for intermediate files  

## Pipeline Steps

Refer to ../excluded_families_testset/README.md.

## Intermediate Files

- Unique miRNA family count file: `unique_family_counts.tsv`

## Output Files

- Final files: 
    - `{primary_input_file}.excluded.tsv`
    - `{primary_input_file}.remaining.tsv`
- Log file: `postprocess_1_exclude_mirna_families.log`