#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=mirbench_test
#SBATCH --ntasks=1

# Initialize variables
OUT_DIR=""

# Parse command-line options
while getopts ":o:" opt; do
  case $opt in
    o) OUT_DIR="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2
      exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2
      exit 1 ;;
  esac
done

# Check if OUT_DIR is provided
if [ -z "$OUT_DIR" ]; then
    echo "Error: Output directory not specified."
    echo "Usage: sbatch run_benchmark_all.sh -o <output_directory>"
    exit 1
fi

python /data/biogemt/miRNA_Binding_Site_Prediction/mirbench/benchmark_all.py --out_dir $OUT_DIR

