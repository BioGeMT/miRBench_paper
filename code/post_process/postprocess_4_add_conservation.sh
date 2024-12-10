#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=pp_4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30

# parse command-line arguments
while getopts i:o:p:c: flag; do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
        p) phyloP_file=${OPTARG};;
        c) phastCons_file=${OPTARG};;
    esac
done

# Check if required arguments are provided
if [ -d "$input_dir" ] || [ -d "$output_dir" ] || [ -f "$phyloP_file" ] || [ -f "$phastCons_file" ]; then
    echo "Usage: $0 -i input_dir -o output_dir -p phyloP_file -c phastCons_file"
    exit 1
fi

# define log file in the output directory
log_file="$output_dir/postprocess_4_conservation.log"

# redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# check for permission issues with log file
if [ $? -ne 0 ]; then
    echo "Error setting up log file. Check your permissions."
    exit 1
fi

# define paths to the directories where the scripts are located
conservation_dir="../conservation"

# define constants for suffixes with extensions
CONSERVATION_SUFFIX=".conservation"

# process each .tsv file in the input directory
for input_file in "$input_dir"/*.tsv; do

    # Add validated conservation scores to all datasets
    base_name=$(basename "$input_file" .tsv)
    input_file="$input_dir/${base_name}.tsv"
    conservation_file="$output_dir/${base_name}${CONSERVATION_SUFFIX}.tsv"
    
    echo "Adding validated conservation scores to $input_file..."
    python3 "$conservation_dir/add_conservation_scores.py" --ifile "$input_file" --phyloP "$phyloP_file" --phastCons "$phastCons_file" --ofile "$conservation_file"
    if [ $? -ne 0 ]; then
        echo "Error in adding conservation scores. Check your script and input file."
        exit 1
    fi
    echo "Conservation scores added. Output saved to $conservation_file"

done 

# Done
echo "Conservation scores successfully added to all datasets in $input_dir."
