import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import pandas as pd


def parse_cutoff(path: Path) -> float:
    match = re.search(r"clusters_cutoff_(.+)\.csv$", path.name)
    if not match:
        raise ValueError(f"Could not parse cutoff from filename: {path.name}")
    return float(match.group(1).replace("p", "."))


def summarize_cluster_file(path: Path) -> dict:
    df = pd.read_csv(path)

    required_cols = {"Gene_ID", "Cluster_ID"}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"{path.name} is missing required columns: {sorted(missing_cols)}")

    cluster_sizes = df["Cluster_ID"].value_counts()

    n_sequences = int(len(df))
    n_clusters = int(cluster_sizes.shape[0])
    n_singletons = int((cluster_sizes == 1).sum())
    largest_cluster_size = int(cluster_sizes.max())

    return {
        "cutoff": parse_cutoff(path),
        "file": path.name,
        "n_sequences": n_sequences,
        "n_clusters": n_clusters,
        "mean_cluster_size": float(cluster_sizes.mean()),
        "median_cluster_size": float(cluster_sizes.median()),
        "min_cluster_size": int(cluster_sizes.min()),
        "max_cluster_size": largest_cluster_size,
        "n_singletons": n_singletons,
        "singleton_fraction": n_singletons / n_clusters if n_clusters else 0.0,
        "largest_cluster_fraction": largest_cluster_size / n_sequences if n_sequences else 0.0,
    }


def make_plots(summary_df: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_specs = [
        ("n_clusters", "Number of Clusters", "cutoff_vs_n_clusters"),
        ("singleton_fraction", "Singleton Fraction", "cutoff_vs_singleton_fraction"),
        ("max_cluster_size", "Max Cluster Size", "cutoff_vs_max_cluster_size"),
    ]

    for column, y_label, stem in plot_specs:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(summary_df["cutoff"], summary_df[column], marker="o", linewidth=2)
        ax.set_xlabel("Cutoff")
        ax.set_ylabel(y_label)
        ax.set_title(f"Cutoff vs {y_label}")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        png_path = output_dir / f"{stem}.png"
        fig.savefig(png_path, dpi=300)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize cutoff sweep cluster outputs.")
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing clusters_cutoff_*.csv files",
    )
    parser.add_argument(
        "--output_file",
        default=None,
        help="Path to write the summary table",
    )
    parser.add_argument(
        "--plot_dir",
        default=None,
        help="Directory to write summary plots",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    if input_dir.name != "genes_cutoff_sweep":
        raise ValueError(
            f"Expected --input_dir to end with 'genes_cutoff_sweep', got: {input_dir}"
        )
    output_root = input_dir.parent

    if args.output_file is None:
        args.output_file = str(output_root / "genes_cutoff_sweep_summary.tsv")
    if args.plot_dir is None:
        args.plot_dir = str(output_root / "summarize_cutoff_plots")

    cluster_files = sorted(input_dir.glob("clusters_cutoff_*.csv"), key=parse_cutoff)

    if not cluster_files:
        raise FileNotFoundError(f"No clusters_cutoff_*.csv files found in {input_dir}")

    summary_df = pd.DataFrame(summarize_cluster_file(path) for path in cluster_files)
    summary_df = summary_df.sort_values("cutoff").reset_index(drop=True)
    summary_df.to_csv(args.output_file, sep="\t", index=False)
    make_plots(summary_df, Path(args.plot_dir))

    print(f"Saved summary to {args.output_file}")
    print(f"Saved plots to {args.plot_dir}")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
