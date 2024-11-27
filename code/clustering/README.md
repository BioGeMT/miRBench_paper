# Gene Sequence Clustering Pipeline

Pipeline for clustering gene sequences and mapping results back to source data.

## Steps

### 1. Generate FASTA (`gene_fasta.py`)

**Input:** TSV file with 'gene' column containing sequences
**Output:** FASTA file with sequences (>Seq_N headers)

```bash
python gene_fasta.py --input data.tsv --output sequences.fasta
```

### 2. Cluster Sequences (`clustering.R`)

**Input:** FASTA file from step 1
**Output:** CSV with sequence IDs and cluster assignments
**Parameters:** Clustering cutoff = 0.1


### 3. Map Clusters (`map_gene_clusters.py`)

**Input:**
- FASTA file from step 1
- Cluster CSV from step 2
- Original TSV file

**Output:** TSV file with original data plus ClusterID column

```bash
python map_gene_clusters.py \
  --fasta_file sequences.fasta \
  --cluster_csv clusters_cutoff_0.1.csv \
  --gene_tsv data.tsv \
  --output_tsv final_output.tsv
```

## Requirements
- Python: pandas, biopython
- R: BiocManager, Biostrings, DECIPHER