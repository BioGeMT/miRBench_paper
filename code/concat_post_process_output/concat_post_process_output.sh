#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=HD_postprocess_concat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Parse options
while getopts "i:o:" opt; do
  case $opt in
    i) input_dir=$OPTARG ;;
    o) output_dir=$OPTARG ;;
    \?) echo "Usage: $0 -i <input_directory> -o <output_directory>"
        exit 1 ;;
  esac
done

# Check if both input and output directories are provided
if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
  echo "Usage: $0 -i <input_directory> -o <output_directory>"
  exit 1
fi

# Define log file in the output directory
log_file="$output_dir/HD_postprocess_concat.log"

# Redirect all output to the log file and also to the terminal
exec > >(tee -a "$log_file") 2>&1

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Define arrays for the types and sizes
types=("train_data" "test_data")
sizes=("1" "10" "100")

# Loop through each type (train_data, test_data)
for type in "${types[@]}"; do
  # Loop through each size (1, 10, 100)
  for size in "${sizes[@]}"; do
    output_file="${output_dir}/AGO2_eCLIP_Manakov2022_${size}_${type}_dataset.tsv"
    
    # Find all matching files
    files=(${input_dir}/*_${type}_with_negatives_${size}.tsv)
    
    # Write header from the first file to the output file
    head -n 1 "${files[0]}" > "$output_file"
    
    # Append content from all files, skipping the header line
    for file in "${files[@]}"; do
      tail -n +2 "$file" >> "$output_file"
    done
    
    echo "Created ${output_file}"
  done
done
