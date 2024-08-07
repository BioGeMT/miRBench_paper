#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=mirbench_test
#SBATCH --ntasks=1

# Initialize variables
OUT_DIR=""
DWNLD_DIR=""

# Parse command-line options
while getopts ":o:d:" opt; do
  case $opt in
    o) OUT_DIR="$OPTARG" ;;
    d) DWNLD_DIR="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2
      exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2
      exit 1 ;;
  esac
done

# Check if OUT_DIR is provided
if [ -z "$OUT_DIR" ]; then
    echo "Error: Output directory not specified. This is where the results will be saved."
    echo "Usage: sbatch run_benchmark_all.sh -o <output_directory>"
    exit 1
fi

# Check if DWNLD_DIR is provided
if [ -z "$DWNLD_DIR" ]; then
    echo "Error: Download directory not specified. This is where the datasets will be downloaded."
    echo "Usage: sbatch run_benchmark_all.sh -o <output_directory> -d <download_directory>"
    exit 1
fi

python benchmark_all.py --out_dir $OUT_DIR --download_dir $DWNLD_DIR

