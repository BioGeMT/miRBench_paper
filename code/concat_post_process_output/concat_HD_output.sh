#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=concat_HD_output
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Parse command line arguments
output_dir="output"  # Default value for output_dir
output_prefix="concat"  # Default value for output_prefix

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_dir) input_dir="$2"; shift ;;
        --output_dir) output_dir="$2"; shift ;;
        --output_prefix) output_prefix="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check if required arguments are provided
if [ -z "$input_dir" ]; then
  echo "Usage: $0 --input_dir <input_directory> [--output_dir <output_directory> (default: output)] [--output_prefix <output_prefix> (default: concat)]"
  exit 1
fi

# Check if output directory exists, else create it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# Full path for output file
output_file="${output_prefix}.unified_length_all_types_unique_high_confidence.tsv"
output_path="$output_dir/$output_file"

# Full path for log file
log_file="$output_dir/concat_HD_output.log"

# Redirect all stdout and stderr to the log file and also to the terminal
exec > >(tee -a "$log_file") 2>&1

# Initialize a flag to check if we've written the header
header_written=false

# Loop through all matching files
for file in "$input_dir"/*.unified_length_all_types_unique_high_confidence.tsv; do
  if [ -f "$file" ] || { [ -L "$file" ] && [ -f "$(readlink -f "$file")" ]; }; then  
    if [ "$header_written" = false ]; then
        # Write the header from the first file
        head -n1 "$file" > "$output_path"
        header_written=true
    fi
    # Append the content (excluding header) to the output file
    tail -n+2 "$file" >> "$output_path"
  else
    echo "Skipping $file: Not a regular file or valid symlink."
  fi
done

echo -e "Concatenation of all *.unified_length_all_types_unique_high_confidence.tsv files in $input_dir directory complete. \nOutput saved to $output_path."

