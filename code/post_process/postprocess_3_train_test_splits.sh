#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=pp_3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30

# parse command-line arguments
while getopts i:o: flag; do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

# Check if required arguments are provided
if [ -d "$input_dir" ] || [ -d "$output_dir" ]; then
    echo "Usage: $0 -i input_dir -o output_dir"
    exit 1
fi

# define log file in the output directory
log_file="$output_dir/postprocess_3_train_test_split.log"

# redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# check for permission issues with log file
if [ $? -ne 0 ]; then
    echo "Error setting up log file. Check your permissions."
    exit 1
fi

# define constants for suffixes with extensions
TRAIN_SUFFIX=".train"
TEST_SUFFIX=".test"

# process each .tsv file in the input directory
for input_file in "$input_dir"/*.negatives.tsv; do

    # Split Train/Test based on the test (6th) column
    echo "Splitting data into train and test sets based on the test column..."

    base_name=$(basename "$input_file" .tsv)

    train_file="$output_dir/${base_name}${TRAIN_SUFFIX}.tsv"
    test_file="$output_dir/${base_name}${TEST_SUFFIX}.tsv"

    awk -F'\t' 'NR==1{header=$0; print header > "'"$train_file"'"; print header > "'"$test_file"'"} NR>1{if($6=="False"){print > "'"$train_file"'"} else {print > "'"$test_file"'"}}' "$input_file"
    if [ $? -ne 0 ]; then
        echo "Error in splitting data. Check your input file."
        exit 1
    fi

    echo "Data split for $input_file completed. Train set saved to $train_file, test set saved to $test_file"

done

# Done
echo "Train and test split pipeline completed successfully."