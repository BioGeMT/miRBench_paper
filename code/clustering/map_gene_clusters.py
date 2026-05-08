import argparse
import pandas as pd

def main(cluster_csv, dataset_tsv, lookup_tsv, output_tsv):
    
    clusters_df = pd.read_csv(cluster_csv)
    gene_df = pd.read_csv(dataset_tsv, sep="\t")
    lookup_df = pd.read_csv(lookup_tsv, sep="\t")

    if "Gene_ID" not in clusters_df.columns:
        raise ValueError("Cluster CSV must contain a 'Gene_ID' column.")

    merged_df = gene_df.merge(
        lookup_df,
        on="gene",
        how="left",
        validate="many_to_one",
    )

    if merged_df["gene_id"].isna().any():
        missing_count = int(merged_df["gene_id"].isna().sum())
        raise ValueError(f"Failed to map gene IDs for {missing_count} rows.")

    merged_df = merged_df.merge(
        clusters_df[["Gene_ID", "Cluster_ID"]],
        left_on="gene_id",
        right_on="Gene_ID",
        how="left",
        validate="many_to_one",
    )

    if merged_df["Cluster_ID"].isna().any():
        missing_count = int(merged_df["Cluster_ID"].isna().sum())
        raise ValueError(f"Failed to map cluster IDs for {missing_count} rows.")

    merged_df = merged_df.drop(columns=["gene_id", "Gene_ID"])
    merged_df = merged_df.rename(columns={"Cluster_ID": "gene_cluster_ID"})
    
    merged_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Results saved to {output_tsv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map cluster IDs to sequences and merge with gene data.")
    parser.add_argument("--cluster_csv", required=True, help="CSV file containing cluster information")
    parser.add_argument("--dataset_tsv", required=True, help="Input dataset TSV file")
    parser.add_argument("--lookup_tsv", required=True, help="TSV file mapping gene sequences to integer gene IDs")
    parser.add_argument("--output_tsv", required=True, help="Output TSV file path")
    args = parser.parse_args()
    main(args.cluster_csv, args.dataset_tsv, args.lookup_tsv, args.output_tsv)
