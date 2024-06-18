#!/bin/bash

# define input, intermediate, and output directories
input_dir="input_data"
intermediate_dir="intermediate_data"
output_dir="output_data"

# create directories if they don't exist
mkdir -p $input_dir $intermediate_dir $output_dir

# define log file in the output directory
log_file="$output_dir/pipeline.log"

# redirect all output to the log file
exec > >(tee -a $log_file) 2>&1

# define input files
input_file="$input_dir/data.tsv"
mature_file="$input_dir/mature.fa"

# step 1: filtering
echo "running filtering step..."
python3 filtering.py --ifile $input_file --ofile $intermediate_dir/filtered_data.tsv
echo "filtering completed. output saved to $intermediate_dir/filtered_data.tsv"

# step 2: family assignment
echo "running family assignment step..."
python3 family_assign.py --ifile $intermediate_dir/filtered_data.tsv --mature $mature_file --ofile $intermediate_dir/family_assigned_data.tsv
echo "family assignment completed. output saved to $intermediate_dir/family_assigned_data.tsv"

# step 3: split train/test based on the test column
echo "splitting data into train and test sets based on the test column..."
train_file="$intermediate_dir/train_data.tsv"
test_file="$intermediate_dir/test_data.tsv"
awk -F'\t' 'NR==1{header=$0; print header > "'$train_file'"; print header > "'$test_file'"} NR>1{if($5=="False"){print > "'$train_file'"} else {print > "'$test_file'"}}' $intermediate_dir/family_assigned_data.tsv
echo "data split completed. train set saved to $train_file, test set saved to $test_file"

# step 4: make negatives for train and test sets with different ratios
for ratio in 1 10 100; do
    echo "generating negative samples with ratio $ratio for train set..."
    python3 make_neg_sets.py --ifile $train_file --ofile $output_dir/train_data_with_negatives_$ratio.tsv --neg_ratio $ratio
    echo "negative samples with ratio $ratio for train set saved to $output_dir/train_data_with_negatives_$ratio.tsv"

    echo "generating negative samples with ratio $ratio for test set..."
    python3 make_neg_sets.py --ifile $test_file --ofile $output_dir/test_data_with_negatives_$ratio.tsv --neg_ratio $ratio
    echo "negative samples with ratio $ratio for test set saved to $output_dir/test_data_with_negatives_$ratio.tsv"
done

# done
echo "post-processing pipeline completed successfully."
