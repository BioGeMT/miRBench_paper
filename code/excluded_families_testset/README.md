# Creation of dataset with excluded unique miRNA families of one dataset

These scripts find all unique miRNA families of one dataset (compared against 2 others) and create two new datasets originating from the original one excluding or including these families. 

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
* Does not account for families with the value "0" or "unknown"
* Output includes both the family name and occurrence count

### 2. Family Filter (`dataset_split_based_on_unique_families.py`)

Creates two datasets by splitting the input based on specific ncRNA families.

**Input:**
* Original dataset (TSV format)
* File containing families to filter by. The noncodingRNA_fam values of this file will be selected and removed from the original file into a new dataset (see excluded dataset output). The remaining original dataset after this removal will be saved (see remaining dataset output)
* Both files must contain a 'noncodingRNA_fam' column

**Output:**
* Two TSV files:
  * One containing rows with the unique noncodingRNA_fam values that were removed from the original file (excluded dataset)
  * One containing all the remaining rows after the above removal (remaining dataset)
* Prints summary statistics showing row counts for both outputs

**Usage:**
```bash
python dataset_split_based_on_unique_families.py \
  --input original_dataset.tsv \
  --families allowed_families.tsv \
  --excluded_dataset excluded_dataset.tsv \
  --remaining_dataset leftout_dataset.tsv
```

## Requirements
* Python 3.x
* pandas
* argparse

## File Format Requirements

All input TSV files must contain a column named 'noncodingRNA_fam' that contains the family annotations.
