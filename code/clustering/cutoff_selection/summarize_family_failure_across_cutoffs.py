import argparse
import re
from pathlib import Path

import pandas as pd


def parse_cutoff(path: Path) -> float:
    match = re.search(r"clusters_cutoff_(.+)\.mapped\.sorted\.tsv$", path.name)
    if not match:
        raise ValueError(f"Could not parse cutoff from filename: {path.name}")
    return float(match.group(1).replace("p", "."))


def summarize_family_for_file(mapped_sorted_path: Path, family: str) -> dict:
    df = pd.read_csv(mapped_sorted_path, sep="\t")
    fam = df[df["noncodingRNA_fam"] == family].copy()

    all_clusters = set(df["gene_cluster_ID"].astype(int).unique().tolist())
    fam_clusters = set(fam["gene_cluster_ID"].astype(int).unique().tolist())
    eligible_clusters = all_clusters - fam_clusters

    requested_negs = int(len(fam))
    eligible_negative_clusters = int(len(eligible_clusters))
    slack = eligible_negative_clusters - requested_negs

    return {
        "cutoff": parse_cutoff(mapped_sorted_path),
        "file": mapped_sorted_path.name,
        "family": family,
        "n_rows": requested_negs,
        "n_unique_noncodingRNA_name": int(fam["noncodingRNA_name"].nunique()),
        "n_unique_noncodingRNA": int(fam["noncodingRNA"].nunique()),
        "n_unique_gene_cluster_ID": int(fam["gene_cluster_ID"].nunique()),
        "n_unique_clusters_total": int(len(all_clusters)),
        "requested_negs": requested_negs,
        "eligible_negative_clusters": eligible_negative_clusters,
        "slack": int(slack),
        "fails": bool(slack < 0),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize family-level negative-sampling feasibility across cutoff mapped TSVs.",
    )
    parser.add_argument(
        "--mapped_dir",
        required=True,
        help="Directory containing clusters_cutoff_*.mapped.sorted.tsv files.",
    )
    parser.add_argument(
        "--family",
        default="mir-17",
        help="Family label in noncodingRNA_fam to summarize.",
    )
    parser.add_argument(
        "--output_path",
        default=None,
        help="Path to write summary TSV. Defaults to <family>_family_summary.tsv in current working directory.",
    )
    args = parser.parse_args()

    mapped_dir = Path(args.mapped_dir)
    if args.output_path is None:
        output_path = Path.cwd() / f"{args.family}_family_summary.tsv"
    else:
        output_path = Path(args.output_path)

    mapped_files = sorted(
        mapped_dir.glob("clusters_cutoff_*.mapped.sorted.tsv"),
        key=parse_cutoff,
    )
    if not mapped_files:
        raise FileNotFoundError(f"No mapped sorted files found in {mapped_dir}")

    rows = [summarize_family_for_file(path, args.family) for path in mapped_files]
    summary_df = pd.DataFrame(rows).sort_values("cutoff").reset_index(drop=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(output_path, sep="\t", index=False)

    print(f"Saved summary to {output_path}")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()