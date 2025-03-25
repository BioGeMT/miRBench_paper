# Gene Cluster Analysis

This script analyses the cluster IDs assigned to the intermediate file with `.gene_clusters_added` suffix from ../post_process/postprocess_2_make_negatives.sh, performing three main tasks:

1. Plots Cluster Size Distribution

2. Computes Biggest Cluster ID, Size, and Percentage of total genes in this cluster 

3. Computes Singleton Cluster Percentage of total genes in singleton clusters


## Requirements

- Python 3
- pandas
- matplotlib

## Usage

Run the script from the command line using:

```bash
python analyse_gene_clusters.py --input <input_file.tsv> --plot_output <plot_output.png> --metrics_output <metrics_output.txt>
```
