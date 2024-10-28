# Note on datasets used for getting plots

## Hejret2023 dataset
For *Hejret2023* dataset, you can find the input file [here](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/HybriDetector/refs/heads/main/ML/Datasets/AGO2_CLASH_Hejret2023_full_dataset.tsv). 
You need to filter for 'miRNA' only in the 'noncodingRNA_type' column and do deduplication on 'seq.g' + 'noncodingRNA_seq' while keeping the first line occurrence. You then need to change 'miRNA_fam' column name to 'noncodingRNA_fam', and 'noncodingRNA_seq' column name to 'noncodingRNA', and run miRBench_paper/code/family_assign/family_assign.py to add miRNA family where value is 0. 

## Klimentova2022 dataset
For *Klimentova2022* dataset, merge [this raw HybriDetector output](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/miRBind/refs/heads/main/Datasets/AGO2_eCLIP_Klimentova22_full_dataset.tsv)
filtered for 'miRNA' only in the 'noncodingRNA_type' column with
[the file prepared for model testing](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/miRBind/refs/heads/main/Datasets/AGO2_eCLIP_Klimentova2022_1_test.tsv) filtered for positives only ('label' == 1). 
Do inner merge based on gene sequence ('seq.g' and 'gene') and miRNA sequence ('noncodingRNA_seq' and 'noncodingRNA') columns. You then need to change 'miRNA_fam' column name to 'noncodingRNA_fam', and 'noncodingRNA_seq' column name to 'noncodingRNA', and run miRBench_paper/code/family_assign/family_assign.py to add miRNA family where value is 0. 

## Manakov2022 dataset
For *Manakov2022* dataset, use the intermediate output from the post-process pipeline from the family_assign step (before the 'sorting' and 'make_negatives' steps) with the `.unified_length_all_types_unique_high_confidence_family_assigned_data.tsv` suffix. Code for plotting uses different column names than in the Manakov dataset. Change the column names in the Manakov dataset as follows:
gene / noncodingRNA / noncodingRNA_fam to seq.g / seq.m / miRNA_fam, respectively. 
