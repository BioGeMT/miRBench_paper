# postprocess_6_genomic_region_annotation_README

## Usage

Can be submitted to SLURM workload manager with `sbatch` or:
```bash
postprocess_6_genomic_region_annotation.sh -i input_dir -o output_dir
```

### Parameters

`-i`: Input directory containing .tsv files with chr, start, end, strand columns  
`-o`: Output directory (step1/ and step2/ subdirectories will be created automatically)

## Pipeline Steps

1. **annotate**: Overlaps intervals against Ensembl release 90 transcripts and produces per-nt region matrix
2. **summarize-sites**: Selects best transcript per interval using CLASH policy and generates region summaries

Refer to https://github.com/BioGeMT/genomic_region_annotator for full documentation.

## Output Files

- Step 1: `output_dir/step1/{basename}.annotated_transcripts.tsv`, `{basename}.annotated_matrix.tsv`, `{basename}.annotated_step1_stats.tsv`
- Step 2: `output_dir/step2/{basename}.annotated_site_summary.tsv` — Main output with all original input columns plus all annotation columns
- Step 2 (Supporting): `{basename}.annotated_final.tsv`, `{basename}.annotated_step2_stats.tsv`
- Log file: `postprocess_6_genomic_region_annotation.log`

## Added Columns (from Step 2 summarize-sites)

The main output TSV (`*_site_summary.tsv`) preserves all original input columns and adds the following annotation columns:

**Transcript Selection:**
- `selected_transcript_id`: Ensembl transcript ID chosen by CLASH policy
- `selected_gene_id`: Ensembl gene ID of selected transcript
- `selected_gene_name`: HGNC gene name of selected transcript

**Region Classification (Dominant):**
- `dominant_region_selected`: Single best genomic region for the selected transcript (UTR3, CDS, UTR5, EXON_OTHER, INTRON, INTERGENIC)
- `dominant_region_union`: Single best region across all passing transcripts (union evidence)

**Region Classification (Multi-region):**
- `regions_present_selected`: Pipe-separated list of all regions present in selected transcript (e.g., "CDS|UTR3")
- `regions_present_union`: Pipe-separated list of all regions present in union evidence

**Base Pair Counts (selected transcript):**
- `bp_utr3_selected`, `bp_cds_selected`, `bp_utr5_selected`, `bp_exon_other_selected`, `bp_intron_selected`, `bp_intergenic_selected`

**Base Pair Counts (union evidence):**
- `bp_utr3_union`, `bp_cds_union`, `bp_utr5_union`, `bp_exon_other_union`, `bp_intron_union`, `bp_intergenic_union`

**Ambiguity and Evidence:**
- `ambiguous_union_vs_selected`: Flag (0/1) indicating whether union evidence differs from selected transcript dominant region
- `n_passing_transcripts`: Number of transcripts with evidence for this interval

## Main Final Output

The primary output file is `{basename}.annotated_site_summary.tsv` in `output_dir/step2/`. This file contains all original input columns (id, gene, noncodingRNA, test, label, etc.) plus the complete set of annotation results:
- Selected transcript information (transcript_id, gene_id, gene_name)
- Dominant region classifications (both selected transcript and union evidence)
- Base pair overlap counts for each genomic region (UTR3, CDS, UTR5, EXON_OTHER, INTRON, INTERGENIC)
- Ambiguity flags and transcript evidence counts

Use this file for downstream analysis and model training. The stats files (`*_step1_stats.tsv`, `*_step2_stats.tsv`) provide summaries of the annotation process and region composition patterns.
