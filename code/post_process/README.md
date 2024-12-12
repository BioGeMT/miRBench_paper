
# Post-Processing Pipelines

This series of pipelines is designed to process as input, the HybriDetector `*.unified_length_all_types_unique_high_confidence.tsv` output files. 


It is intended to be used on the following datasets:
- https://github.com/ML-Bioinfo-CEITEC/HybriDetector/blob/main/ML/Datasets/AGO2_CLASH_Hejret2023_full_dataset.tsv 
- https://github.com/ML-Bioinfo-CEITEC/miRBind/blob/main/Datasets/AGO2_eCLIP_Klimentova22_full_dataset.tsv
- the concatenated output of HD when processing selected samples from the Manakov data (to be uploaded somewhere still)

Note that for the Hejret and Klimentova datasets above, the `miRNA_fam` column must be renamed to `noncodingRNA_fam` after downloading, prior to any processing, for consistency of all column names.

The scope of the pipeline series is to filter the files for miRNA data, deduplicate gene-miRNA sequence pairs, create a left-out test set with miRNA families unique only to this set, construct the negative class in an unbiased manner, split the datasets into training and testing, and finally add conservation score to the gene sequences. 


The series is composed of 6 pipelines (listed below) and are intended to be run in the defined order as the output of one feeds the next. Refer to the worflow diagram. 

1. postprocess_0_filter_and_deduplicate
2. postprocess_1_exclude_mirna_families
3. postprocess_2_make_negatives
4. postprocess_3_train_test_splits
5. postprocess_4_drop_test_col
6. postprocess_4_add_conservation

## Requirements
- Python 3
- Run `conda env create --file=post_process.yml`, then `conda activate postprocess`
- Manually download the `hg38.phyloP100way.bw` and `hg38.phastCons100way.bw` files from:
   - https://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/ 
   - https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/
- Ensure the necessary Python scripts are located in the specified relative paths:
  - `../filtering/filtering.py`
  - `../excluded_families_testset/unique_family_counter.py`
  - `../excluded_families_testset/dataset_split_based_on_unique_families.py`
  - `../clustering/gene_fasta.py`
  - `../clustering/clustering.R `
  - `../clustering/map_gene_clusters.py`
  - `../make_neg_sets/make_neg_sets.py`
  - `../conservation/add_conservation_scores.py`

## Usage

Each pipeline of the series must be run separately. Refer to the corresponding README files. 
