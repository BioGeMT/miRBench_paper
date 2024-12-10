import argparse
import pandas as pd

def main(cluster_csv, dataset_tsv, output_tsv):
    
    clusters_df = pd.read_csv(cluster_csv)
    gene_df = pd.read_csv(dataset_tsv, sep="\t")
    
    gene_df["gene_cluster_ID"] = clusters_df["Cluster_ID"]
    
    gene_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Results saved to {output_tsv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map cluster IDs to sequences and merge with gene data.")
    parser.add_argument("--cluster_csv", required=True, help="CSV file containing cluster information")
    parser.add_argument("--dataset_tsv", required=True, help="Input dataset TSV file")
    parser.add_argument("--output_tsv", required=True, help="Output TSV file path")
    args = parser.parse_args()
    main(args.cluster_csv, args.dataset_tsv, args.output_tsv)
