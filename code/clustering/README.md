# Gene Sequence Clustering Pipeline

A pipeline for clustering gene sequences and mapping cluster assignments back to source data.

## Overview

This pipeline consists of three main steps:
1. Convert TSV data to FASTA format
2. Perform sequence clustering using DECIPHER
3. Map cluster assignments back to the original dataset

## Pipeline Steps

### 1. Generate FASTA (`gene_fasta.py`)

Converts a TSV file containing gene sequences into FASTA format for clustering on unique target sequences only.

**Input:** 
- TSV file with a 'gene' column containing sequences

**Output:** 
- FASTA file containing unique gene sequences with integer `gene_id` headers
- TSV lookup file mapping `gene_id` back to the original `gene` sequence

**Usage:**
```bash
python gene_fasta.py --input data.tsv --output sequences.fasta
```

### 2. Cluster Sequences (`clustering.R`)

Performs sequence clustering using the DECIPHER package.

**Input:** 
- FASTA file from step 1

**Output:** 
- CSV file containing integer gene IDs and their cluster assignments

**Parameters:**
- Clustering cutoff: 0.1
- Number of processors: 8

**Usage:**
```bash
Rscript clustering.R sequences.fasta clusters.csv
```

### 3. Map Clusters (`map_gene_clusters.py`)

Maps cluster assignments back to the original dataset.

**Input:**
- Cluster CSV file from step 2
- Original dataset TSV file
- Lookup TSV file from step 1

**Output:**
- TSV file containing original data plus a new 'gene_cluster_ID' column

Cluster IDs are mapped back in two steps: first from `gene` to `gene_id` using the lookup TSV, then from `gene_id` to `Cluster_ID` using the clustering output. This allows clustering to be performed on unique targets while assignments are propagated to all matching rows in the original dataset.

**Usage:**
```bash
python map_gene_clusters.py \
  --cluster_csv clusters.csv \
  --dataset_tsv data.tsv \
  --lookup_tsv sequences.gene_id_lookup.tsv \
  --output_tsv final_output.tsv

```

## Requirements

### Python Dependencies
- pandas
- argparse

### R Dependencies
- BiocManager
- Biostrings (Bioconductor)
- DECIPHER (Bioconductor)
