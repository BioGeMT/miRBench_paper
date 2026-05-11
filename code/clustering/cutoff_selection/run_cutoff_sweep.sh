#!/bin/bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
default_clustering_script="$(cd "$script_dir/.." && pwd)/clustering.R"

usage() {
    cat <<EOF
Usage: $0 -i input_fasta [-o output_dir] [-c cutoff_list] [-p processors] [-r clustering_script]

Required:
    -i  Input FASTA file

Optional:
  -o  Output directory for cluster CSV files
      (default: ./outputs/<dataset_stem>/genes_cutoff_sweep, inferred from -i)
  -c  Comma-separated cutoff list (default: 0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3)
  -p  Number of processors to pass to clustering.R (default: 8)
  -r  Path to clustering.R
      (default: ../clustering.R relative to this script)
EOF
}

clustering_script="$default_clustering_script"
input_fasta=""
output_dir=""
cutoff_csv="0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3"
processors=8

while getopts "i:o:c:p:r:h" flag; do
    case "${flag}" in
        i) input_fasta=${OPTARG} ;;
        o) output_dir=${OPTARG} ;;
        c) cutoff_csv=${OPTARG} ;;
        p) processors=${OPTARG} ;;
        r) clustering_script=${OPTARG} ;;
        h)
            usage
            exit 0
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done

if [ -z "$input_fasta" ]; then
    echo "Input FASTA is required (-i)."
    usage
    exit 1
fi

if [ ! -f "$input_fasta" ]; then
    echo "Input FASTA not found: $input_fasta"
    exit 1
fi

if [ -z "$output_dir" ]; then
    input_fasta_name="$(basename "$input_fasta")"
    dataset_stem="${input_fasta_name%.genes.fasta}"
    dataset_stem="${dataset_stem%.fasta}"
    output_dir="$script_dir/outputs/$dataset_stem/genes_cutoff_sweep"
fi

if [ ! -f "$clustering_script" ]; then
    echo "Clustering script not found: $clustering_script"
    exit 1
fi

mkdir -p "$output_dir"

IFS=',' read -r -a cutoffs <<< "$cutoff_csv"

for cutoff in "${cutoffs[@]}"; do
    sanitized_cutoff=${cutoff//./p}
    output_file="$output_dir/clusters_cutoff_${sanitized_cutoff}.csv"

    echo "Running clustering for cutoff=$cutoff"
    Rscript "$clustering_script" "$input_fasta" "$output_file" "$cutoff" "$processors"
    echo "Saved $output_file"
done

echo "Cutoff sweep completed."
