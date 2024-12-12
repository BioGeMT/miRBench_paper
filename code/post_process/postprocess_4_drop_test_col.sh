#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=pp_4
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
if [ ! -d "$input_dir" ] || [ ! -d "$output_dir" ]; then
    echo "Usage: $0 -i input_dir -o output_dir"
    exit 1
fi

# define log file in the output directory
log_file="$output_dir/postprocess_4_drop_test_col.log"

# redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# check for permission issues with log file
if [ $? -ne 0 ]; then
    echo "Error setting up log file. Check your permissions."
    exit 1
fi

# define constants for suffixes with extensions
DROP_TEST_COL_SUFFIX=".drop_test_col"

# process each .tsv file in the input directory
for input_file in "$input_dir"/*.tsv; do

    # Removing the 'test' (6th) column from the input file
    echo "Removing the 'test' (6th) column from $input_file"

    base_name=$(basename "$input_file" .tsv)
    drop_test_col_file="$output_dir/${base_name}${DROP_TEST_COL_SUFFIX}.tsv"

    # Use awk to remove the 6th column without causing column shifts
    awk -F'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) if(i!=6) printf "%s%s", $i, (i==NF?"\n":OFS)}' "$input_file" > "$drop_test_col_file"

    if [ $? -ne 0 ]; then
        echo "Error in removing the 'test' (6th) column."
        exit 1
    fi   

    echo "'test' column removed successfully. Output saved to $drop_test_col_file"

done

# Done
echo "Removing 'test' column for all files in input dir completed successfully."