#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=HD_postprocess_on_concat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

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
default_ratios=(1 10 100)
neg_ratios=( "${neg_ratios[@]:-"${default_ratios[@]}"}" )
min_edit_distance=${min_edit_distance:-3}

# define paths to the directories where the scripts are located
filtering_dir="../filtering"
family_assign_dir="../family_assign"
make_neg_sets_dir="../make_neg_sets"

# define directories for output and intermediate files
output_dir="${output_dir:-$(pwd)/output}"
intermediate_dir="${intermediate_dir:-$(pwd)/intermediate}"
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
DEDUPLICATED_SUFFIX="_deduplicated_data.tsv"
FAMILY_ASSIGNED_SUFFIX="_family_assigned_data.tsv"
SORTED_FAMILY_ASSIGNED_SUFFIX="_family_assigned_data_sorted.tsv"
TRAIN_SUFFIX="_train_"
TEST_SUFFIX="_test_"
NEG_SUFFIX="_with_negatives_"

# process each .tsv file in the input directory
for input_file in "$input_dir"/*unified_length_all_types_unique_high_confidence.tsv; do
    base_name=$(basename "$input_file" .tsv)
    filtered_file="$intermediate_dir/${base_name}${FILTERED_SUFFIX}"
    deduplicated_file="$intermediate_dir/${base_name}${DEDUPLICATED_SUFFIX}"
    family_assigned_file="$intermediate_dir/${base_name}${FAMILY_ASSIGNED_SUFFIX}"
    family_assigned_file_sorted="$intermediate_dir/${base_name}${SORTED_FAMILY_ASSIGNED_SUFFIX}"

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

    # Step 3: Family Assignment
    echo "Running family assignment step on $deduplicated_file..."
    python3 "$family_assign_dir/family_assign.py" --ifile "$deduplicated_file" --mature "$mature_file" --ofile "$family_assigned_file"
    if [ $? -ne 0 ]; then
        echo "Error in family assignment step. Check your script and input file."
        exit 1
    fi
    echo "Family assignment completed. Output saved to $family_assigned_file"

    # Step 4: Sort the family assigned file based on the first column in preparation for negative sample generation
    echo "Sorting the family assigned file based on the first column..."
    (head -n 1 "$family_assigned_file" && tail -n +2 "$family_assigned_file" | sort -k 1) > "${family_assigned_file_sorted}"
    echo "Family assigned file sorted. Output saved to $family_assigned_file_sorted"

    # Step 5: Make negatives with different ratios
    for ratio in ${neg_ratios[@]}; do
        echo "Generating negative samples with ratio $ratio and min edit distance $min_edit_distance..."

        neg_output="$intermediate_dir/${base_name}${NEG_SUFFIX}${ratio}.tsv"

        echo "Running make_neg_sets.py for $family_assigned_file_sorted with ratio $ratio..."
        python3 "$make_neg_sets_dir/make_neg_sets.py" --ifile "$family_assigned_file_sorted" --ofile "$neg_output" --neg_ratio "$ratio" --min_required_edit_distance "$min_edit_distance"
        if [ $? -ne 0 ]; then
            echo "Error in generating negative samples. Check your script and input file."
            exit 1
        fi
        echo "File with negative samples for ratio $ratio saved to $neg_output"
    done
    echo "Negative samples generation completed for $input_file"

    # Step 6: Split Train/Test based on the test column
    echo "Splitting data into train and test sets based on the test column..."
    for ratio in ${neg_ratios[@]}; do
        neg_file="$intermediate_dir/${base_name}${NEG_SUFFIX}${ratio}.tsv"
        train_file="$output_dir/${base_name}${TRAIN_SUFFIX}${ratio}.tsv"
        test_file="$output_dir/${base_name}${TEST_SUFFIX}${ratio}.tsv"
        awk -F'\t' 'NR==1{header=$0; print header > "'"$train_file"'"; print header > "'"$test_file"'"} NR>1{if($5=="False"){print > "'"$train_file"'"} else {print > "'"$test_file"'"}}' "$neg_file"
        if [ $? -ne 0 ]; then
            echo "Error in splitting data. Check your input file."
            exit 1
        fi
        echo "Data split completed. Train set saved to $train_file, test set saved to $test_file"
    done

    # Step 7: Remove the fifth column from the train and test files
    echo "Removing the fifth column from the train and test sets..."
    for ratio in "${neg_ratios[@]}"; do
        for suffix in "$TRAIN_SUFFIX" "$TEST_SUFFIX"; do
            file="$output_dir/${base_name}${suffix}${ratio}.tsv"
            # Use awk to remove the fifth column without causing column shifts
            awk -F'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) if(i!=5) printf "%s%s", $i, (i==NF?"\n":OFS)}' "$file" > "${file}_tmp" && mv "${file}_tmp" "$file"
        done
    done
    if [ $? -ne 0 ]; then
        echo "Error in removing the fifth column."
        exit 1
    fi
    
    echo "Fifth column removed from train set and test set."
    
    done

# Done
echo "Post-processing pipeline completed successfully."
