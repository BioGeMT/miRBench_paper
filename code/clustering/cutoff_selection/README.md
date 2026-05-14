# Cutoff Selection Analysis

This directory contains the full cutoff-sweep analysis pipeline used to evaluate clustering behavior and downstream negative-sampling feasibility.

The workflow has 3 stages:
1. Run clustering across multiple cutoffs.
2. Summarize cluster structure metrics per cutoff.
3. Analyze downstream block feasibility per cutoff.

## Files

- `run_cutoff_sweep.sh`
  - Runs `clustering.R` for each cutoff and writes `clusters_cutoff_*.csv`.
- `summarize_cutoff_sweep.py`
  - Summarizes clustering metrics and generates summary plots.
- `analyze_downstream_utility.py`
  - Maps clusters to dataset rows, builds negative-sampling blocks, and summarizes failure/slack metrics.

## Prerequisites

From `miRBench_paper/code`:
- `clustering/clustering.R`
- `clustering/map_gene_clusters.py`
- `make_neg_sets/make_neg_sets.py`
- `sort_by_column/sort_tsv.sh`

Runtime requirements:
- Python 3 with `pandas` and `matplotlib`
- R with `Rscript` available for `clustering.R`

## Expected Input Artifacts

- Gene FASTA used for clustering (for step 1)
  - Produced by the `gene_fasta` script.
- Gene lookup TSV mapping unique gene sequence to integer `Gene_ID` (for step 3)
  - Produced by the `gene_fasta` script.
- Positive dataset TSV in post-process format (for step 3)
  - For Hejret and Klimentova: output of post-process 0.
  - For Manakov excluded and remaining: output of post-process 1.

## Step 1: Run Cutoff Sweep

Run from `miRBench_paper/code/clustering/cutoff_selection`:

```bash
bash run_cutoff_sweep.sh \
  -i <path/to/dataset.genes.fasta> \
  -p 8
```

Optional flags:
- `-o` output directory (default: `./outputs/<dataset_stem>/genes_cutoff_sweep`)
- `-c` comma-separated cutoffs (default: `0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3`)
- `-r` path to clustering script (default: `../clustering.R`)

Output:
- `outputs/<dataset_stem>/genes_cutoff_sweep/clusters_cutoff_*.csv`

## Step 2: Summarize Cluster Structure

```bash
python summarize_cutoff_sweep.py \
  --input_dir outputs/<dataset_stem>/genes_cutoff_sweep
```

Notes:
- `--input_dir` must end with `genes_cutoff_sweep`.

Outputs (default):
- `outputs/<dataset_stem>/genes_cutoff_sweep_summary.tsv`
- `outputs/<dataset_stem>/summarize_cutoff_plots/`
  - `cutoff_vs_n_clusters.png`
  - `cutoff_vs_singleton_fraction.png`
  - `cutoff_vs_max_cluster_size.png`

## Step 3: Analyze Downstream Utility

```bash
python analyze_downstream_utility.py \
  --dataset_tsv <path/to/positive_dataset.tsv> \
  --cluster_dir outputs/<dataset_stem>/genes_cutoff_sweep \
  --lookup_tsv <path/to/lookup.tsv>
```

Default outputs:
- `outputs/<dataset_stem>/genes_cutoff_sweep_downstream_utility.tsv`
- `outputs/<dataset_stem>/downstream_utility_plots/`
- `outputs/<dataset_stem>/intermediate_mapped_tsvs/`
  - persistent per-cutoff mapped files used for block analysis:
  - `clusters_cutoff_*.mapped.sorted.tsv`

Important behavior:
- Mapped unsorted TSVs are temporary only.
- Mapped sorted TSVs are persisted to `intermediate_mapped_tsvs`.
- The downstream summary includes:
  - failure counts and fractions
  - `failing_block_labels`
  - slack statistics
  - available negative-cluster statistics

## Output Layout (Default)

```text
outputs/<dataset_stem>/
  genes_cutoff_sweep/
    clusters_cutoff_0p01.csv
    ...
  genes_cutoff_sweep_summary.tsv
  summarize_cutoff_plots/
    cutoff_vs_n_clusters.png
    cutoff_vs_singleton_fraction.png
    cutoff_vs_max_cluster_size.png
  genes_cutoff_sweep_downstream_utility.tsv
  downstream_utility_plots/
    cutoff_vs_failing_blocks.png
    cutoff_vs_median_slack.png
    cutoff_vs_min_slack.png
    cutoff_vs_median_available_negative_clusters.png
    cutoff_vs_min_available_negative_clusters.png
  intermediate_mapped_tsvs/
    clusters_cutoff_0p01.mapped.sorted.tsv
    ...
```

## Optional: Targeted Failing-Family Analysis

This was used to report the per-cutoff failure pattern for mir-17 in Manakov remaining. Use this helper to report the per-cutoff failure pattern when one family is failing across cutoffs.

Script:
- `summarize_family_failure_across_cutoffs.py`

Example:

```bash
python summarize_family_failure_across_cutoffs.py \
  --mapped_dir outputs/<dataset_stem>/intermediate_mapped_tsvs \
  --family mir-17 \
  --output_path outputs/<dataset_stem>/mir-17_family_summary.tsv
```

Inputs:
- `--mapped_dir` should point to the mapped sorted files produced in step 3.

Outputs:
- family summary TSV with requested negatives, eligible negative clusters, slack, and fail flag per cutoff

