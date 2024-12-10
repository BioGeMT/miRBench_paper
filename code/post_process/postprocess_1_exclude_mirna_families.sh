#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=pp_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30

# parse command-line arguments
while getopts i:j:k:o:n: flag; do
    case "${flag}" in
        i) input_file=${OPTARG};;
        j) input_file_2=${OPTARG};;
        k) input_file_3=${OPTARG};;
        o) output_dir=${OPTARG};;
        n) intermediate_dir=${OPTARG};;
    esac
done

# check if required argument is provided
if [ -f "$input_file" ] || [ -f "$input_file_2" ] || [ -f "$input_file_3" ] || [ -d "$output_dir" ] || [ -d "$intermediate_dir" ]; then
    echo "Usage: $0 -i input_file -j input_file_2 -k input_file_3 -o output_dir -n intermediate_dir"
    exit 1
fi

# define log file in the output directory
log_file="$output_dir/postprocess_1_exclude_mirna_families.log"

# redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# check for permission issues with log file
if [ $? -ne 0 ]; then
    echo "Error setting up log file. Check your permissions."
    exit 1
fi

# define paths to the directories where the scripts are located
exclude_families_dir="../excluded_families_testset"

# define constants for suffixes with extensions
EXCLUDED_SUFFIX=".excluded"
REMAINING_SUFFIX=".remaining"

# process the input file to get the base name
base_name=$(basename "$input_file" .tsv)

# define paths to the intermediate and output files
counts_file="$intermediate_dir/unique_family_counts.tsv"
excluded_file="$output_dir/${base_name}${EXCLUDED_SUFFIX}.tsv"
remaining_file="$output_dir/${base_name}${REMAINING_SUFFIX}.tsv"

# Step 1: Identifying and counting miRNA families that are unique in one dataset relative to another two datasets
echo "Running unique family counts step on $input_file..."
python3 "$exclude_families_dir/unique_family_counter.py" --unique "$input_file" --file2 "$input_file_2" --file3 "$input_file_3" --output "$counts_file"
if [ $? -ne 0 ]; then
    echo "Error in running unique family counts step. Check your script and input file."
    exit 1
fi
echo "Unique family counting completed. Output saved to $counts_file"

# Step 2: Filtering dataset based on miRNA families that are unique in one dataset relative to another two datasets
echo "Running filtering unique families step on $input_file..."
python3 "$exclude_families_dir/dataset_split_based_on_unique_families.py" --input "$input_file" --families "$counts_file" --excluded_dataset "$excluded_file" --remaining_dataset "$remaining_file"
if [ $? -ne 0 ]; then
    echo "Error in filtering unique families step. Check your script and input file."
    exit 1
fi
echo "Filtering unique families step completed. Outputs saved to $excluded_file and $remaining_file"

# Done
echo "Identifying and filtering unique miRNA families in one dataset relative to another two datasets pipeline completed successfully."