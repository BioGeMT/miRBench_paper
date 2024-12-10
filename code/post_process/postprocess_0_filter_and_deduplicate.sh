#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=pp_0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30

# parse command-line arguments
while getopts i:o:n: flag; do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
        n) intermediate_dir=${OPTARG};;
    esac
done

# check if required argument is provided
if [ -d "$input_dir" ] || [ -d "$output_dir" ] || [ -d "$intermediate_dir" ]; then
    echo "Usage: $0 -i input_dir -o output_dir -n intermediate_dir"
    exit 1
fi

# define log file in the output directory
log_file="$output_dir/postprocess_0_filter_and_deduplicate.log"

# redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# check for permission issues with log file
if [ $? -ne 0 ]; then
    echo "Error setting up log file. Check your permissions."
    exit 1
fi

# define paths to the directories where the scripts are located
filtering_dir="../filtering"

# define constants for suffixes with extensions
FILTERED_SUFFIX=".filtered"
DEDUPLICATED_SUFFIX=".deduplicated"

# process each .tsv file in the input directory
for input_file in "$input_dir"/*.tsv; do
    base_name=$(basename "$input_file" .tsv)
    filtered_file="$intermediate_dir/${base_name}${FILTERED_SUFFIX}.tsv"
    deduplicated_file="$output_dir/${base_name}${FILTERED_SUFFIX}${DEDUPLICATED_SUFFIX}.tsv"

    # Step 1: Filtering
    echo "Running filtering step on $input_file..."
    python3 "$filtering_dir/filtering.py" --ifile "$input_file" --ofile "$filtered_file"
    if [ $? -ne 0 ]; then
        echo "Error in filtering step. Check your script and input file."
        exit 1
    fi
    echo "Filtering completed. Output saved to $filtered_file"

    # Step 2: Deduplication
    echo "Running deduplication step on $filtered_file..."
    # deduplicate based on combination of first two columns
    awk -F'\t' 'NR==1{print $0} NR>1{if(!seen[$1$2]++){print}}' "$filtered_file" > "$deduplicated_file"
    if [ $? -ne 0 ]; then
        echo "Error in deduplication step. Check your input file."
        exit 1
    fi
    echo "Deduplication completed. Output saved to $deduplicated_file"
done

# Done
echo "Filter and deduplication pipeline completed successfully."