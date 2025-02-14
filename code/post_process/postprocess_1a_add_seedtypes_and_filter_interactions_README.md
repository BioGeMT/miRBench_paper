# postprocess_1a_add_seed_types_and_filter_interactions_README

## Usage

Can be submitted to SLURM workload manager with `sbatch` or:
```bash
postprocess_1a_add_seed_types_and_filter_interactions.sh -i input_dir -o output_dir -n intermediate_dir
```
### Parameters

`-i`: Input directory containing .tsv files  
`-o`: Output directory for final files  
`-n`: Intermediate directory for intermediate files  

## Pipeline Steps

1. **Adding seed types**: Refer to ../filter_interactions/add_seed_types_README.md
2. **Filtering interactions**: Refer to ../filter_interactions/filter_interactions_README.md

## Intermediate Files

- Input file with added seed types (Seed6mer, Seed6merBulgeOrMismatch): `{input_filename}.seed_types.tsv`

## Output Files

- Input file filtered for canonical seed: `{input_filename}.canonical6mer.tsv`
- Input file filtered for non-canonical seed: `{input_filename}.noncanonical6mer.tsv`
- Input file filtered for non-seed interactions: `{input_filename}.nonseed.tsv`
- Log file: `postprocess_1a_add_seedtypes_and_filter_interactions.log`