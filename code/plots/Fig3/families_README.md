
## Description

The `families.py` script analyzes miRNA family prevalence and seed type distribution across multiple datasets. It processes the input data, calculates family counts and percentages, and generates both a summary TSV file and a plot visualizing the seed type distributions for the top miRNA families.

## Usage

Run the script from the command line with the following syntax:

```
python families.py dataset1.tsv dataset2.tsv dataset3.tsv -o output_plot.png -t output_summary.tsv
```

### Arguments:
- `dataset1.tsv`, `dataset2.tsv`, `dataset3.tsv`: Input TSV files containing miRNA data.
- `-o`, `--output`: (Optional) The name of the output PNG file for the plot. Default is 'seed_prevalence_plot.png'.
- `-t`, `--tsv_output`: (Optional) The name of the output TSV file for the summary results. Default is 'miRNA_family_analysis.tsv'.

## Input File Format

The input TSV files should have a header row and include columns named 'seq.g' (guide sequence), 'seq.m' (miRNA sequence), and 'miRNA_fam' (miRNA family).

For *Hejret2023* dataset, you can find the input file [here](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/HybriDetector/refs/heads/main/ML/Datasets/AGO2_CLASH_Hejret2023_full_dataset.tsv). 
You need to filter for 'miRNA' only in the 'noncodingRNA_type' column and do deduplication on 'seq.g' + 'noncodingRNA_seq' while keeping the first line occurrence.

For *Klimentova2022* dataset, merge [this raw HybriDetector output](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/miRBind/refs/heads/main/Datasets/AGO2_eCLIP_Klimentova22_full_dataset.tsv)
filtered for 'miRNA' only in the 'noncodingRNA_type' column with
[the file prepared for model testing](https://raw.githubusercontent.com/ML-Bioinfo-CEITEC/miRBind/refs/heads/main/Datasets/AGO2_eCLIP_Klimentova2022_1_test.tsv) filtered for positives only ('label' == 1). 
Do inner merge based on gene sequence ('seq.g' and 'gene') and miRNA sequence ('noncodingRNA_seq' and 'noncodingRNA') columns.

For *Manakov2022* use product of the pipeline right before the 'make negative' step.


## Output

The script generates two output files:
1. A TSV file containing counts and percentages of miRNA families across all datasets.
2. A PNG file with a grouped bar plot showing the distribution of seed types for the top 10 miRNA families.

## Dependencies

- pandas
- matplotlib
- argparse
- Bio (Biopython)

Make sure to have `seed_utils.py` in the same directory as this script.
