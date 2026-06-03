# postprocess_2_make_negatives_README

## Usage

Can be submitted to SLURM workload manager with `sbatch` or:
```bash
postprocess_2_make_negatives.sh -i input_dir -o output_dir -n intermediate_dir
```
### Parameters

`-i`: Input directory containing .tsv files  
`-o`: Output directory for final files  
`-n`: Intermediate directory for intermediate files  

## Pipeline Steps

1. **Clustering**: Refer to ../clustering/README.md
2. **Generating negatives**: Refer to ../make_neg_sets/README.md

## Intermediate Files

- Clustering files: 
    - FASTA file: `{input_file}.fasta`
    - Clustering output file: `{input_file}.gene_clusters.csv`
    - Input file with cluster ID: `{input_file}.gene_clusters_added.tsv`

## Output Files

- Final file: `{input_filename}.gene_clusters_added.negatives.tsv`
- Log file: `postprocess_2_make_negatives.log`