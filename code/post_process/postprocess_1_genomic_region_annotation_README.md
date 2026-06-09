# postprocess_1_genomic_region_annotation_README

## Usage

Can be submitted to SLURM workload manager with `sbatch` or:
```bash
postprocess_1_genomic_region_annotation.sh -i input_dir -o output_dir
```

### Parameters

`-i`: Input directory containing .tsv files with chr, start, end, strand columns  
`-o`: Output directory (`step1/` and `step2/` subdirectories will be created automatically)

## Pipeline Steps

1. **annotate**: Overlaps intervals against Ensembl release 90 transcripts and produces per-nt region matrix
2. **summarize-sites**: Selects best transcript per interval using CLASH policy and generates region summaries
3. **filter annotated summaries**: Filters the Step 2 `*.annotated_site_summary.tsv` output to the downstream columns needed for later processing and model input

Refer to https://github.com/BioGeMT/genomic_region_annotator for full documentation.

## Output Files

- Step 1: `output_dir/step1/{basename}.annotated_transcripts.tsv`, `{basename}.annotated_matrix.tsv`, `{basename}.annotated_step1_stats.tsv`
- Step 2: `output_dir/step2/{basename}.annotated_site_summary.tsv` — Full genomic region annotation output with all original input columns plus annotation results
- Step 2 (Supporting): `{basename}.annotated_final.tsv`, `{basename}.annotated_step2_stats.tsv`
- Step 3: `output_dir/{basename}.annotated.filtered.tsv` — Filtered downstream output with selected columns only
- Log file: `postprocess_6_genomic_region_annotation.log`

## Filtered Final Output

The final downstream output file is `{basename}.annotated.filtered.tsv` in `output_dir/`. It is created from the Step 2 file `{basename}.annotated_site_summary.tsv` and keeps only the following columns in this order:

- `gene`
- `noncodingRNA`
- `noncodingRNA_name`
- `noncodingRNA_fam`
- `feature`
- `test`
- `label`
- `chr`
- `start`
- `end`
- `strand`
- `Nunique`
- `dominant_region`
- `regions_present`
- `read_start_in_sel_tx_1based`
- `read_end_in_sel_tx_1based`

The final four columns are renamed from the Step 2 site summary columns:

- `dominant_region_selected` -> `dominant_region`
- `regions_present_selected` -> `regions_present`
- `selected_read_start_in_tx_1based` -> `read_start_in_sel_tx_1based`
- `selected_read_end_in_tx_1based` -> `read_end_in_sel_tx_1based`

Use `{basename}.annotated.filtered.tsv` for downstream analysis and model training. Use the full Step 2 site summary file when you need the complete annotation output, and the stats files (`*_step1_stats.tsv`, `*_step2_stats.tsv`) for summaries of the annotation process and region composition patterns.
