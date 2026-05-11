#!/bin/bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
default_clustering_script="$(cd "$script_dir/.." && pwd)/clustering.R"
dataset_stem="AGO2_eCLIP_Klimentova22_full_dataset.filtered.deduplicated"
default_input_fasta="$(cd "$script_dir/.." && pwd)/${dataset_stem}.genes.fasta"
default_output_dir="$script_dir/outputs/${dataset_stem}_genes_cutoff_sweep"

usage() {
    cat <<EOF
Usage: $0 -i input_fasta -o output_dir [-c cutoff_list] [-p processors] [-r clustering_script]

Optional:
  -i  Input FASTA file
      (default: ../AGO2_eCLIP_Klimentova22_full_dataset.filtered.deduplicated.genes.fasta relative to this script)
  -o  Output directory for cluster CSV files
      (default: ./outputs/AGO2_eCLIP_Klimentova22_full_dataset.filtered.deduplicated_genes_cutoff_sweep)
  -c  Comma-separated cutoff list (default: 0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3)
  -p  Number of processors to pass to clustering.R (default: 8)
  -r  Path to clustering.R
      (default: ../clustering.R relative to this script)
EOF
}

clustering_script="$default_clustering_script"
input_fasta="$default_input_fasta"
output_dir="$default_output_dir"
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

if [ ! -f "$input_fasta" ]; then
    echo "Input FASTA not found: $input_fasta"
    exit 1
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
