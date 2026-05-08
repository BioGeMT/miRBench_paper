import argparse
import pandas as pd


def require_columns(df, required, name):
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"{name} is missing required columns: {sorted(missing)}")


def main(cluster_csv, dataset_tsv, lookup_tsv, output_tsv):
    clusters = pd.read_csv(cluster_csv)
    dataset = pd.read_csv(dataset_tsv, sep="\t")
    lookup = pd.read_csv(lookup_tsv, sep="\t")

    require_columns(clusters, ["Gene_ID", "Cluster_ID"], "Cluster CSV")
    require_columns(dataset, ["gene"], "Dataset TSV")
    require_columns(lookup, ["gene", "gene_id"], "Lookup TSV")

    if lookup["gene"].duplicated().any():
        dupes = lookup.loc[lookup["gene"].duplicated(), "gene"].head().tolist()
        raise ValueError(f"Lookup TSV has duplicate gene entries, e.g. {dupes}")

    if clusters["Gene_ID"].duplicated().any():
        dupes = clusters.loc[clusters["Gene_ID"].duplicated(), "Gene_ID"].head().tolist()
        raise ValueError(f"Cluster CSV has duplicate Gene_ID entries, e.g. {dupes}")

    gene_to_id = lookup.set_index("gene")["gene_id"]
    id_to_cluster = clusters.set_index("Gene_ID")["Cluster_ID"]

    dataset["gene_id"] = dataset["gene"].map(gene_to_id)

    missing_gene_ids = dataset["gene_id"].isna()
    if missing_gene_ids.any():
        examples = dataset.loc[missing_gene_ids, "gene"].drop_duplicates().head().tolist()
        raise ValueError(
            f"Failed to map gene IDs for {missing_gene_ids.sum()} rows. "
            f"Examples: {examples}"
        )

    dataset["gene_cluster_ID"] = dataset["gene_id"].map(id_to_cluster)

    missing_clusters = dataset["gene_cluster_ID"].isna()
    if missing_clusters.any():
        examples = dataset.loc[missing_clusters, "gene_id"].drop_duplicates().head().tolist()
        raise ValueError(
            f"Failed to map cluster IDs for {missing_clusters.sum()} rows. "
            f"Example gene_ids: {examples}"
        )

    dataset = dataset.drop(columns=["gene_id"])
    dataset.to_csv(output_tsv, sep="\t", index=False)

    print(f"Results saved to {output_tsv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Map gene cluster IDs back onto a dataset TSV."
    )
    parser.add_argument("--cluster_csv", required=True)
    parser.add_argument("--dataset_tsv", required=True)
    parser.add_argument("--lookup_tsv", required=True)
    parser.add_argument("--output_tsv", required=True)
    args = parser.parse_args()

    main(args.cluster_csv, args.dataset_tsv, args.lookup_tsv, args.output_tsv)