# Creation of dataset with excluded unique miRNA families

These scripts find all unique miRNA families of one dataset (compared against 2 others) and create a new dataset including only those unique rows. 

## Scripts

### 1. Unique Family Counter (`unique_family_counter.py`)

Identifies and counts families that are unique to one dataset when compared against two others.

**Input:**
* Primary TSV file to analyze
* Two additional TSV files to compare against
* All files must contain a 'noncodingRNA_fam' column

**Output:**
* TSV file with unique families and their counts
* Prints summary statistics to console

**Usage:**
```bash
python unique_family_counter.py \
  --unique primary_data.tsv \
  --file2 compare_data1.tsv \
  --file3 compare_data2.tsv \
  --output unique_families.tsv
```

**Notes:**
* Automatically excludes 'unknown' and '0' family annotations
* Output includes both the family name and occurrence count

### 2. Family Filter (`exclude_fam.py`)

Filters a dataset to keep only entries matching specific ncRNA families.

**Input:**
* Original dataset (TSV format)
* File containing allowed families (from unique_family_counter.py output)
* Both files must contain a 'noncodingRNA_fam' column

**Output:**
* Filtered TSV file containing only rows with matching families

**Usage:**
```bash
python exclude_fam.py \
  --input original_dataset.tsv \
  --families allowed_families.tsv \
  --output filtered_dataset.tsv
```

## Requirements
* Python 3.x
* pandas
* argparse

## File Format Requirements

All input TSV files must contain a column named 'noncodingRNA_fam' that contains the family annotations.
