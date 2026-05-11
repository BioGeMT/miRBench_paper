import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))
sys.path.append(str(Path(__file__).resolve().parents[2] / "make_neg_sets"))
from map_gene_clusters import map_clusters_to_dataset
from make_neg_sets import get_negative_sampling_block_label, yield_negative_sampling_sub_blocks


def parse_cutoff(path: Path) -> float:
    match = re.search(r"clusters_cutoff_(.+)\.csv$", path.name)
    if not match:
        raise ValueError(f"Could not parse cutoff from filename: {path.name}")
    return float(match.group(1).replace("p", "."))


def iter_negative_sampling_blocks(df: pd.DataFrame):
    for _, fam_block in df.groupby("noncodingRNA_fam", sort=False):
        for block in yield_negative_sampling_sub_blocks(fam_block):
            yield block.reset_index(drop=True)


def load_cluster_assignments(cluster_path: Path) -> pd.DataFrame:
    clusters_df = pd.read_csv(cluster_path)

    required_cols = {"Gene_ID", "Cluster_ID"}
    missing_cols = required_cols - set(clusters_df.columns)
    if missing_cols:
        raise ValueError(
            f"{cluster_path.name} is missing required columns: {sorted(missing_cols)}.")

    if clusters_df["Gene_ID"].duplicated().any():
        dup_count = int(clusters_df["Gene_ID"].duplicated().sum())
        raise ValueError(f"{cluster_path.name} has {dup_count} duplicate Gene_ID values.")

    return clusters_df[["Gene_ID", "Cluster_ID"]]


def summarize_blocks(mapped_df: pd.DataFrame) -> dict:
    all_clusters = set(mapped_df["gene_cluster_ID"].unique().tolist())
    block_rows = []

    for block in iter_negative_sampling_blocks(mapped_df):
        block_clusters = set(block["gene_cluster_ID"].unique().tolist())
        num_neg = int(block.shape[0])
        positive_cluster_count = len(block_clusters)
        allowed_clusters = all_clusters - block_clusters
        available_negative_clusters = len(allowed_clusters)
        slack = available_negative_clusters - num_neg

        block_rows.append(
            {
                "block_label": get_negative_sampling_block_label(block),
                "block_size": num_neg,
                "positive_cluster_count": positive_cluster_count,
                "allowed_cluster_count": len(allowed_clusters),
                "available_negative_clusters": available_negative_clusters,
                "slack": slack,
                "fails": slack < 0,
            }
        )

    blocks_df = pd.DataFrame(block_rows)
    failing_blocks = blocks_df[blocks_df["fails"]]

    return {
        "n_rows": int(mapped_df.shape[0]),
        "n_unique_clusters_total": int(len(all_clusters)),
        "n_blocks": int(blocks_df.shape[0]),
        "mean_block_size": float(blocks_df["block_size"].mean()),
        "median_block_size": float(blocks_df["block_size"].median()),
        "mean_positive_clusters_per_block": float(blocks_df["positive_cluster_count"].mean()),
        "median_positive_clusters_per_block": float(blocks_df["positive_cluster_count"].median()),
        "mean_allowed_clusters_per_block": float(blocks_df["allowed_cluster_count"].mean()),
        "median_allowed_clusters_per_block": float(blocks_df["allowed_cluster_count"].median()),
        "mean_available_negative_clusters": float(blocks_df["available_negative_clusters"].mean()),
        "median_available_negative_clusters": float(blocks_df["available_negative_clusters"].median()),
        "min_available_negative_clusters": int(blocks_df["available_negative_clusters"].min()),
        "mean_slack": float(blocks_df["slack"].mean()),
        "median_slack": float(blocks_df["slack"].median()),
        "min_slack": int(blocks_df["slack"].min()),
        "n_failing_blocks": int(failing_blocks.shape[0]),
        "failure_fraction": float(failing_blocks.shape[0] / blocks_df.shape[0]),
    }


def make_plots(summary_df: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_specs = [
        ("n_failing_blocks", "Failing Blocks", "cutoff_vs_failing_blocks"),
        ("median_slack", "Median Slack", "cutoff_vs_median_slack"),
        ("min_slack", "Minimum Slack", "cutoff_vs_min_slack"),
        ("median_available_negative_clusters", "Median Available Negative Clusters", "cutoff_vs_median_available_negative_clusters"),
    ]

    for column, y_label, stem in plot_specs:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(summary_df["cutoff"], summary_df[column], marker="o", linewidth=2)
        ax.set_xlabel("Cutoff")
        ax.set_ylabel(y_label)
        ax.set_title(f"Cutoff vs {y_label}")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(output_dir / f"{stem}.png", dpi=300)
        plt.close(fig)


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    dataset_stem = "AGO2_eCLIP_Klimentova22_full_dataset.filtered.deduplicated"

    parser = argparse.ArgumentParser(description="Analyze downstream negative-sampling utility across clustering cutoffs.")
    parser.add_argument(
        "--dataset_tsv",
        default=str(base_dir.parent / "AGO2_eCLIP_Klimentova22_full_dataset.filtered.deduplicated.tsv"),
        help="Positive dataset TSV in post-process format.",
    )
    parser.add_argument(
        "--cluster_dir",
        default=str(base_dir / "outputs" / f"{dataset_stem}_genes_cutoff_sweep"),
        help="Directory containing clusters_cutoff_*.csv files.",
    )
    parser.add_argument(
        "--lookup_tsv",
        default=str(base_dir.parent / f"{dataset_stem}.genes.gene_id_lookup.tsv"),
        help="TSV file mapping unique gene sequences to integer gene IDs.",
    )
    parser.add_argument(
        "--output_file",
        default=str(base_dir / "outputs" / f"{dataset_stem}_genes_cutoff_sweep_downstream_utility.tsv"),
        help="Path to write the downstream utility summary TSV.",
    )
    parser.add_argument(
        "--plot_dir",
        default=str(base_dir / "outputs" / f"{dataset_stem}_genes_cutoff_sweep" / "downstream_utility_plots"),
        help="Directory to write downstream utility plots.",
    )
    args = parser.parse_args()

    dataset_df = pd.read_csv(args.dataset_tsv, sep="\t")
    lookup_df = pd.read_csv(args.lookup_tsv, sep="\t")
    cluster_paths = sorted(Path(args.cluster_dir).glob("clusters_cutoff_*.csv"), key=parse_cutoff)

    if not cluster_paths:
        raise FileNotFoundError(f"No clusters_cutoff_*.csv files found in {args.cluster_dir}")

    summary_rows = []
    for cluster_path in cluster_paths:
        clusters_df = load_cluster_assignments(cluster_path)
        mapped_df = map_clusters_to_dataset(dataset_df, lookup_df, clusters_df)
        block_summary = summarize_blocks(mapped_df)
        block_summary["cutoff"] = parse_cutoff(cluster_path)
        block_summary["file"] = cluster_path.name
        summary_rows.append(block_summary)

    summary_df = pd.DataFrame(summary_rows).sort_values("cutoff").reset_index(drop=True)
    summary_df.to_csv(args.output_file, sep="\t", index=False)
    make_plots(summary_df, Path(args.plot_dir))

    print(f"Saved downstream utility summary to {args.output_file}")
    print(f"Saved downstream utility plots to {args.plot_dir}")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
