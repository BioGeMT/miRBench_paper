# Post-Processing Pipelines

This series of pipelines is designed to process as input the HybriDetector `*.unified_length_all_types_unique_high_confidence.tsv` output files.

It is intended to be used on the following datasets:
- https://github.com/ML-Bioinfo-CEITEC/HybriDetector/blob/main/ML/Datasets/AGO2_CLASH_Hejret2023_full_dataset.tsv
- https://github.com/ML-Bioinfo-CEITEC/miRBind/blob/main/Datasets/AGO2_eCLIP_Klimentova22_full_dataset.tsv
- https://zenodo.org/records/14501607/files/AGO2_eCLIP_Manakov2022_full_dataset.tsv.gz

Note that for the Hejret and Klimentova datasets above, the `miRNA_fam` column must be renamed to `noncodingRNA_fam` after downloading, prior to any processing, for consistency of all column names.

The scope of the pipeline series is to filter the files for miRNA data, annotate all intervals with genomic region information (UTR, CDS, intron, etc.) using Ensembl transcripts, filter the annotation output to the downstream columns used by later steps, deduplicate gene-miRNA sequence pairs, create a left-out test set with miRNA families unique only to this set, construct the negative class in an unbiased manner, split the datasets into training and testing, and finally add conservation score to the gene sequences.

The series is composed of 7 pipelines (listed below) and are intended to be run in the defined order as the output of one feeds the next. Refer to the workflow diagram.

1. postprocess_0_filter_and_deduplicate
2. postprocess_1_genomic_region_annotation
3. postprocess_2_exclude_mirna_families
4. postprocess_3_make_negatives
5. postprocess_4_train_test_splits
6. postprocess_5_drop_test_col
7. postprocess_6_add_conservation

## Requirements
- Python 3
- Run `conda env create --file=post_process.yml`, then `conda activate postprocess`
- Ensure the required helper scripts and resources are present in the repository:
  - `../filtering/filtering.py`
  - `../genomic_region_annotator_filtering/genomic_region_annotator_filtering.py`
  - `../excluded_families_testset/unique_family_counter.py`
  - `../excluded_families_testset/dataset_split_based_on_unique_families.py`
  - `../clustering/gene_fasta.py`
  - `../clustering/clustering.R`
  - `../clustering/map_gene_clusters.py`
  - `../sort_by_column/sort_tsv.sh`
  - `../make_neg_sets/make_neg_sets.py`
  - `../conservation/add_conservation_scores.py`
- The post-processing shell scripts resolve these helper paths relative to their own location, so they can be run from any working directory.
- For conservation scoring, the required BigWig files are:
  - `hg38.phyloP100way.bw`
  - `hg38.phastCons100way.bw`

## Usage

You can run the post-processing workflow in either of two ways:
- Run each pipeline separately in the documented order. Refer to the corresponding README files.
- Run the wrapper script `run_postprocess_pipeline.sh`.

The wrapper supports two modes:
- `--mode cohort` for the three canonical input datasets in one input directory.
- `--mode single` for one dataset file at a time, including future datasets beyond the canonical three.

Examples:
```bash
bash run_postprocess_pipeline.sh --mode cohort -i input_dir -o output_dir
```

```bash
bash run_postprocess_pipeline.sh --mode single -f input_file.tsv -o output_dir
```

For conservation scoring, the wrapper can use existing files or manage them for you:
- Pass `-p` and `-c` to use existing local BigWig files.
- Omit `-p` and `-c` to let the wrapper use cached files or download them with `wget` into `output_dir/reference_data/conservation/`.

In cohort mode, the wrapper also creates final Zenodo-ready gzipped outputs in `output_dir/zenodo_release/` with the published filenames.

You can also override the cache/download location with:
- `-r <conservation_dir>`
