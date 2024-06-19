#!/bin/bash

# Default directories and parameters
input_dir=""
intermediate_dir="intermediate_data"
output_dir="output_data"
mature_dir="mature_data"
neg_ratios=(1 10 100)
min_edit_distance=3

# Directories for the scripts
filtering_dir="../filtering"
family_assign_dir="../family_assign"
make_neg_sets_dir="../make_neg_sets"

# Parse command-line arguments
while getopts i:o:n:t:r: flag; do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
        n) intermediate_dir=${OPTARG};;
        t) IFS=',' read -r -a neg_ratios <<< "${OPTARG}";;
        r) min_edit_distance=${OPTARG};;
    esac
done

# Check if required argument is provided
if [ -z "$input_dir" ]; then
    echo "Usage: $0 -i input_dir [-o output_dir] [-n intermediate_dir] [-t neg_ratios] [-r min_edit_distance]"
    exit 1
fi

# Create output and intermediate directories if they don't exist
mkdir -p "$intermediate_dir" "$output_dir" "$mature_dir"

# Define log file in the output directory
log_file="$output_dir/pipeline.log"

# Redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# Download mature.fa file if it doesn't exist
mature_file="$mature_dir/mature.fa"
if [ ! -f "$mature_file" ]; then
    echo "Downloading mature.fa file..."
    wget -O "$mature_file" "https://www.mirbase.org/download/mature.fa"
    echo "Download completed. File saved to $mature_file"
fi

# Process each .tsv file in the input directory
for input_file in "$input_dir"/*.tsv; do
    base_name=$(basename "$input_file" .tsv)
    filtered_file="$intermediate_dir/${base_name}_filtered_data.tsv"
    family_assigned_file="$intermediate_dir/${base_name}_family_assigned_data.tsv"
    train_file="$intermediate_dir/${base_name}_train_data.tsv"
    test_file="$intermediate_dir/${base_name}_test_data.tsv"

    # Step 1: Filtering
    echo "Running filtering step on $input_file..."
    python3 "$filtering_dir/filtering.py" --ifile "$input_file" --ofile "$filtered_file"
    echo "Filtering completed. Output saved to $filtered_file"

    # Step 2: Family Assignment
    echo "Running family assignment step on $filtered_file..."
    python3 "$family_assign_dir/family_assign.py" --ifile "$filtered_file" --mature "$mature_file" --ofile "$family_assigned_file"
    echo "Family assignment completed. Output saved to $family_assigned_file"

    # Step 3: Split Train/Test based on the test column
    echo "Splitting data into train and test sets based on the test column..."
    awk -F'\t' 'NR==1{header=$0; print header > "'$train_file'"; print header > "'$test_file'"} NR>1{if($5=="False"){print > "'$train_file'"} else {print > "'$test_file'"}}' "$family_assigned_file"
    echo "Data split completed. Train set saved to $train_file, test set saved to $test_file"

    # Step 4: Make negatives for train and test sets with different ratios
    for ratio in "${neg_ratios[@]}"; do
        echo "Generating negative samples with ratio $ratio and min edit distance $min_edit_distance for train and test sets..."
        python3 "$make_neg_sets_dir/make_neg_sets.py" --ifile "$train_file" --ofile "$output_dir/${base_name}_train_data_with_negatives_${ratio}.tsv" --neg_ratio "$ratio" --min_edit_distance "$min_edit_distance"
        python3 "$make_neg_sets_dir/make_neg_sets.py" --ifile "$test_file" --ofile "$output_dir/${base_name}_test_data_with_negatives_${ratio}.tsv" --neg_ratio "$ratio" --min_edit_distance "$min_edit_distance"
    done
    echo "Negative samples generation completed for $input_file"
done

# Done
echo "Post-processing pipeline completed successfully."
