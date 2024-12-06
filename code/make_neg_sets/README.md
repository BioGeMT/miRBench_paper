# Negative Example Generation Script

## Overview

This script generates **negative examples** by pairing miRNA families with non-overlapping gene clusters. It processes a file containing positive examples and outputs both positive and negative examples to a specified output file.

Negative examples are constructed to ensure they:
- Pair miRNA families with gene clusters not found in the original block and not already picked as a negative for the current block. 
- Match the same number of rows as positive examples, so this script can only produce 1:1 positive to negative class ratio. 

The script processes input data block-wise to optimize memory usage and efficiently handle large datasets.

## Features

- **Block-wise Processing**: Efficiently processes large files by grouping rows based on `miRNA family` and producing negatives per miRNA family.
- **Support for Unknown miRNA Families**: Handles `unknown` miRNA families separately such that negatives are produced per unique miRNA sequence.
- **Negative Example Generation**: Creates negative examples while avoiding cluster overlap with positive or negative examples of the same block.
- **Reproducibility**: Ensures consistent outputs with a fixed random seed.
- **Efficient Output Handling**: Writes processed data directly to the output file to minimize memory usage.

## Prerequisites

- **Python**
- **Required Libraries**:
  - `argparse`
  - `pandas`
  - `random`
  - `time`

## Usage

### Command

python make_neg_sets.py --ifile <input_file> --ofile <output_file>

### Arguments

- `--ifile` (required):  
  Path to the input file containing positive examples. **This file must be sorted by miRNA family.**

- `--ofile` (required):  
  Path to the output file where both positive and negative examples will be saved.

## Input file format

The input file must:

1. Be a tab-separated file (TSV) that is sorted by noncodingRNA_fam to ensure correct block-wise processing. 
2. Contain the following columns (including but not limited to):
    - noncodingRNA_fam: miRNA family identifier.
    - ClusterID: Cluster ID for genes.
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

Negative examples are generated only when there are sufficient non-overlapping clusters. An error is raised otherwise. 
