
# Post-Processing Pipelines

This series of pipelines is designed to process as input, the HybriDetector `*.unified_length_all_types_unique_high_confidence.tsv` output files. 


It is intended to be used on the following dataset:
- https://zenodo.org/records/14501607/files/AGO2_eCLIP_Manakov2022_full_dataset.tsv.gz 

The scope of the pipeline series is to filter the files for miRNA data, deduplicate gene-miRNA sequence pairs, create a left-out test set with miRNA families unique only to this set, annotate seed types and filter canonical, non-canonical, and non-seed interaction types into different files, construct the negative class in an unbiased and class balanced manner, split the datasets into training and testing, and finally add conservation score to the gene sequences. 

The series is composed of 7 pipelines (listed below) and are intended to be run in the defined order as the output of one feeds the next. Refer to the worflow diagram. 

1. postprocess_0_filter_and_deduplicate
2. postprocess_1_exclude_mirna_families
3. postprocess_1a_add_seed_types_and_filter_interactions
4. postprocess_2_make_negatives
5. postprocess_3_train_test_splits
6. postprocess_4_drop_test_col
7. postprocess_5_add_conservation

## Requirements
- Python 3
- Create an environment and make sure the following packages are installed:
  - pandas
  - miRBench (v1.0.0)
  - r-base
  - r-biostrings
  - r-decipher
  - pyBigWig
- Manually download the `hg38.phyloP100way.bw` and `hg38.phastCons100way.bw` files from:
   - https://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/ 
   - https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/
- Ensure the necessary Python scripts are located in the specified relative paths:
  - `../filtering/filtering.py`
  - `../excluded_families_testset/unique_family_counter.py`
  - `../excluded_families_testset/dataset_split_based_on_unique_families.py`
  - `../filter_interactions/add_seed_types.py`
  - `../filter_interactions/filter_interactions.py`
  - `../clustering/gene_fasta.py`
  - `../clustering/clustering.R `
  - `../clustering/map_gene_clusters.py`
  - `../make_neg_sets/make_neg_sets.py`
  - `../conservation/add_conservation_scores.py`

## Usage

Each pipeline of the series must be run separately. Refer to the corresponding README files. 