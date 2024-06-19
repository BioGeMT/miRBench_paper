#!/bin/bash

# parse command-line arguments
while getopts i:o:n:t:r: flag; do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
        n) intermediate_dir=${OPTARG};;
        t) IFS=',' read -r -a neg_ratios <<< "${OPTARG}";;
        r) min_edit_distance=${OPTARG};;
    esac
done

# check if required argument is provided
if [ -z "$input_dir" ]; then
    echo "Usage: $0 -i input_dir [-o output_dir] [-n intermediate_dir] [-t neg_ratios] [-r min_edit_distance]"
    exit 1
fi

# set default values for neg_ratios and min_edit_distance if not specified
neg_ratios=${neg_ratios:-("1" "10" "100")}
min_edit_distance=${min_edit_distance:-3}

# define paths to the directories where the scripts are located
filtering_dir="../filtering"
family_assign_dir="../family_assign"
make_neg_sets_dir="../make_neg_sets"

# define directories for output and intermediate files
output_dir=${output_dir:-$(pwd)/output}
intermediate_dir=${intermediate_dir:-$(pwd)/intermediate}
mature_dir="$intermediate_dir/mature"

# create output and intermediate directories if they don't exist
mkdir -p "$intermediate_dir" "$output_dir" "$mature_dir"

# define log file in the output directory
log_file="$output_dir/pipeline.log"

# redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# check for permission issues with log file
if [ $? -ne 0 ]; then
    echo "Error setting up log file. Check your permissions."
    exit 1
fi

# download mature.fa file if it doesn't exist
mature_file="$mature_dir/mature.fa"
if [ ! -f "$mature_file" ]; then
    echo "Downloading mature.fa file..."
    wget -O "$mature_file" "https://www.mirbase.org/download/mature.fa"
    if [ $? -ne 0 ]; then
        echo "Error downloading mature.fa file. Check your internet connection and permissions."
        exit 1
    fi
    echo "Download completed. File saved to $mature_file"
fi

# define constants for suffixes with extensions
FILTERED_SUFFIX="_filtered_data.tsv"
FAMILY_ASSIGNED_SUFFIX="_family_assigned_data.tsv"
TRAIN_SUFFIX="_train_data.tsv"
TEST_SUFFIX="_test_data.tsv"
NEG_TRAIN_SUFFIX="_train_data_with_negatives_"
NEG_TEST_SUFFIX="_test_data_with_negatives_"

# process each .tsv file in the input directory
for input_file in "$input_dir"/*.tsv; do
    base_name=$(basename "$input_file" .tsv)
    filtered_file="$intermediate_dir/${base_name}${FILTERED_SUFFIX}"
    family_assigned_file="$intermediate_dir/${base_name}${FAMILY_ASSIGNED_SUFFIX}"
    train_file="$intermediate_dir/${base_name}${TRAIN_SUFFIX}"
    test_file="$intermediate_dir/${base_name}${TEST_SUFFIX}"

    # Step 1: Filtering
    echo "Running filtering step on $input_file..."
    python3 "$filtering_dir/filtering.py" --ifile "$input_file" --ofile "$filtered_file"
    if [ $? -ne 0 ]; then
        echo "Error in filtering step. Check your script and input file."
        exit 1
    fi
    echo "Filtering completed. Output saved to $filtered_file"

    # Step 2: Family Assignment
    echo "Running family assignment step on $filtered_file..."
    python3 "$family_assign_dir/family_assign.py" --ifile "$filtered_file" --mature "$mature_file" --ofile "$family_assigned_file"
    if [ $? -ne 0 ]; then
        echo "Error in family assignment step. Check your script and input file."
        exit 1
    fi
    echo "Family assignment completed. Output saved to $family_assigned_file"

    # Step 3: Split Train/Test based on the test column
    echo "Splitting data into train and test sets based on the test column..."
    awk -F'\t' 'NR==1{header=$0; print header > "'$train_file'"; print header > "'$test_file'"} NR>1{if($5=="False"){print > "'$train_file'"} else {print > "'$test_file'"}}' "$family_assigned_file"
    if [ $? -ne 0 ]; then
        echo "Error in splitting data. Check your input file."
        exit 1
    fi
    echo "Data split completed. Train set saved to $train_file, test set saved to $test_file"

    # Step 4: Remove the fifth column from the train and test files
    echo "Removing the fifth column from the train and test sets..."
    awk -F'\t' 'BEGIN{OFS="\t"} { $5=""; sub("\t\t", "\t"); print }' "$train_file" > "${train_file}_tmp" && mv "${train_file}_tmp" "$train_file"
    awk -F'\t' 'BEGIN{OFS="\t"} { $5=""; sub("\t\t", "\t"); print }' "$test_file" > "${test_file}_tmp" && mv "${test_file}_tmp" "$test_file"

    if [ $? -ne 0 ]; then
        echo "Error in removing the fifth column."
        exit 1
    fi

    echo "Fifth column removed from train set and test set."

    # Step 4: Make negatives for train and test sets with different ratios
    for ratio in 1 10 100; do
        echo "Generating negative samples with ratio $ratio and min edit distance $min_edit_distance for train and test sets..."
        neg_train_output="$output_dir/${base_name}${NEG_TRAIN_SUFFIX}${ratio}.tsv"
        neg_test_output="$output_dir/${base_name}${NEG_TEST_SUFFIX}${ratio}.tsv"

        echo "Running make_neg_sets.py for train set with ratio $ratio..."
        python3 "$make_neg_sets_dir/make_neg_sets.py" --ifile "$train_file" --ofile "$neg_train_output" --neg_ratio "$ratio" --min_edit_distance "$min_edit_distance"
        if [ $? -ne 0 ]; then
            echo "Error in generating negative samples for train set. Check your script and input file."
            exit 1
        fi
        echo "Train set with negative samples saved to $neg_train_output"

        echo "Running make_neg_sets.py for test set with ratio $ratio..."
        python3 "$make_neg_sets_dir/make_neg_sets.py" --ifile "$test_file" --ofile "$neg_test_output" --neg_ratio "$ratio" --min_edit_distance "$min_edit_distance"
        if [ $? -ne 0 ]; then
            echo "Error in generating negative samples for test set. Check your script and input file."
            exit 1
        fi
        echo "Test set with negative samples saved to $neg_test_output"
    done
    echo "Negative samples generation completed for $input_file"
done

# Done
echo "Post-processing pipeline completed successfully."
