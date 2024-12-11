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
2. **Sorting by miRNA family**: Sorts TSV file by the 'noncodingRNA_fam' column while preserving header. Automatically finds column position and exits if column not found. Used to prepare data for negative sample generation.
2. **Generating negatives**: Refer to ../make_neg_sets/README.md

## Intermediate Files

- Clustering files: 
    - FASTA file: `{input_file}.fasta`
    - Clustering output file: `{input_file}.gene_clusters.csv`
    - Input file with cluster ID: `{input_file}.gene_clusters_added.tsv`
- Sorted File: `{input_file}.mirfam_sorted.tsv`

## Output Files

- Final file: `{input_filename}.negatives.tsv`
- Log file: `postprocess_2_make_negatives.log`