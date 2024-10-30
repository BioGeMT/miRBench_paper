# Note on datasets used for getting plots

This pipeline prepares the 3 datasets used throughout miRBench_paper/code/plots/. 

## 1. Download datasets

For **Hejret2023** dataset, download [this file](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/HybriDetector/refs/heads/main/ML/Datasets/AGO2_CLASH_Hejret2023_full_dataset.tsv). 

For **Klimentova2022** dataset, download [this raw HybriDetector (HD) output file](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/miRBind/refs/heads/main/Datasets/AGO2_eCLIP_Klimentova22_full_dataset.tsv) and [this file prepared for model testing](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/miRBind/refs/heads/main/Datasets/AGO2_eCLIP_Klimentova2022_1_test.tsv). 

For **Manakov2022** dataset, download [this miRBench_paper post-process pipeline intermediate file](https://zenodo.org/records/13993348/files/concat.unified_length_all_types_unique_high_confidence_with_negatives_1.tsv). You might need to request access from Stephanie Sammut (stephanie.sammut@um.edu.mt). 

## 2. Run `miRBench_paper/code/plots/process_datasets_for_plots.py`

### Requirements

- pandas Python package

### Description

This python script processes the downloaded datasets as follows:

For **Hejret2023** dataset*, it filters the downloaded dataset for miRNAs (column 'noncodingRNA_type' == 'miRNA'), deduplicates the dataset on columns 'seq.g' and 'noncodingRNA_seq', keeping the first line occurrence, and renaming columns (for consistency throughout the miRBench_paper repo) in the order as follows:
- 'seq.g' to 'gene'
- 'noncodingRNA' to 'noncodingRNA_renamed'
- 'noncodingRNA_seq' to 'noncodingRNA'
- 'miRNA_fam' to 'noncodingRNA_fam'
The output is saved to a .tsv file. 

For **Klimentova2022** dataset, it filters the downloaded [raw HD output file](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/miRBind/refs/heads/main/Datasets/AGO2_eCLIP_Klimentova22_full_dataset.tsv) for miRNAs (column 'noncodingRNA_type' == 'miRNA') and renames columns, as for the Hejret2023 dataset above. It then filters the other downloaded [file prepared for model testing](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/miRBind/refs/heads/main/Datasets/AGO2_eCLIP_Klimentova2022_1_test.tsv) filtered for positive examples (column 'label' == 1), and renames the column 'miRNA_fam' to 'noncodingRNA_fam'. The script then does an inner merge of the two processed datasets on columns 'gene' and 'noncodingRNA', and saves the output to a .tsv file. 

For **Manakov2022** dataset, the script filters the downloaded dataset for positive examples (column 'label' == 1), and saves the output to a .tsv file. **This is the file that should be used as input to the scripts under plots/ for the Manakov2022 dataset.**

The script also prints the number of examples in each dataset. This should match the numbers in the miRBench manuscript. 

## 3. Run `miRBench_paper/code//family_assign/family_assign.py` on Hejret2023 and Klimentova2022 .tsv outputs

### Requirements

- pandas Python package
- mature.fa file (Download from miRBase or from miRBench_paper/code/family_assign/family_assign_example_input/mature.fa)

### Description

This step adds the miRNA family where the value is '0'. Refer to miRBench_paper/code/family_assign/family_assign_README.md to run the script (separately) on the output of the Hejret2023 and Klimentova2022 datasets from the previous step. **The output files from this step should be used as input to the scripts under plots/ for the Hejret2023 and Klimentova2022 datasets.**

## Author & Date

Stephanie Sammut
30-October-2024