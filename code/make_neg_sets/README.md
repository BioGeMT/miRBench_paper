# Negative Example Generation Script

## Overview

This script generates **negative examples** by pairing miRNA families with non-overlapping gene clusters paired to other miRNA families. It processes a file containing positive examples and outputs both positive and negative examples to a specified output file.

Negative examples are constructed per miRNA family group to ensure that the abundance of miRNA families in the positive and negative classes are balanced.

Per miRNA family group, we ensure that:
- Gene targets do not belong to the same cluster as the cluster of any positive in that miRNA family group
- Gene targets do not belong to the same cluster from which a gene target has already been selected as a negative in that miRNA family group

This means that the negatives produced for a miRNA family will have the same miRNAs as the positives but are paired with gene targets that are not similar to (do not cluster with) the gene targets of the same miRNAs in the positive class.

Additionally the gene targets from a single cluster are only picked once per miRNA family, ensuring that the miRNA-gene pairs that have the same miRNA family will have gene targets that are not similar to eachother (are not in the same cluster).

This script can only produce 1:1 positive to negative class ratio.

The script processes input data by miRNA-family groups after loading the TSV with pandas, preserving inferred column dtypes across the full dataset and each group.

## Features

- **Block-wise Processing**: Groups rows based on `noncodingRNA_fam` and produces negatives per unique value in the `noncodingRNA_fam` column while preserving pandas-inferred dtypes.
- **Support for Unknown miRNA Families**: Handles `unknown` miRNA families separately such that negatives are produced per unique miRNA sequence.
- **Negative Example Generation**: Creates negative examples with the considerations described in the Overview above.
- **Reproducibility**: Ensures consistent outputs with a fixed random seed.
- **Efficient Output Handling**: Writes processed data directly to the output file to minimize memory usage.

## Prerequisites

- **Python**
- **Required Libraries**:
  - `argparse`
  - `pandas`
  - `time`

## Usage

### Command

```
python make_neg_sets.py --ifile <input_file> --ofile <output_file>
```

### Arguments

- `--ifile` (required):
  Path to the input file containing positive examples.

- `--ofile` (required):
  Path to the output file where both positive and negative examples will be saved.

## Input file format

The input file must:

1. Be a tab-separated file (TSV).
2. Contain the following columns (including but not limited to):
    - noncodingRNA_fam: miRNA family identifier.
    - gene_cluster_ID: Cluster ID for genes.
    - noncodingRNA: miRNA sequence.
    - Additional columns such as gene, feature, chr, etc.

## Output file format

The output file will:

- Contain both positive and negative examples.
- Match the input file's column structure.
- Have a label column with:
    - 1 for positive examples.
    - 0 for negative examples.

## Notes

Negative examples are generated only when there are sufficient candidates to choose from. If a block has fewer eligible negative candidates than positives, the script now warns and randomly downsamples the positives for that block to match the available negative pool. If no eligible negative candidates exist, the block is skipped and the number of excluded positives is reported.
