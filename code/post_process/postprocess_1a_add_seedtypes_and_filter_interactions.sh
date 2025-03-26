#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=pp_1a
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
if [ ! -d "$input_dir" ] || [ ! -d "$output_dir" ] || [ ! -d "$intermediate_dir" ]; then
    echo "Usage: $0 -i input_dir -o output_dir -n intermediate_dir"
    exit 1
fi

# define log file in the output directory
log_file="$output_dir/postprocess_1a_add_seedtypes_and_filter_interactions.log"

# redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# check for permission issues with log file
if [ $? -ne 0 ]; then
    echo "Error setting up log file. Check your permissions."
    exit 1
fi

# define paths to the directories where the scripts are located
add_seeds_and_filtering_dir="../filter_interactions"

# define constants for suffixes with extensions
SEED_TYPES_SUFFIX=".seed_types"
CANONICAL_SUFFIX=".canonical6mer"
NON_CANONICAL_SUFFIX=".noncanonical6mer"
NON_SEED_SUFFIX=".nonseed"

for input_file in "$input_dir"/*.tsv; do
    # process the input file
    base_name=$(basename "$input_file" .tsv)
    seed_types_output="$intermediate_dir/${base_name}${SEED_TYPES_SUFFIX}.tsv"
    canonical_output="$output_dir/${base_name}${CANONICAL_SUFFIX}.tsv"
    noncanonical_output="$output_dir/${base_name}${NON_CANONICAL_SUFFIX}.tsv"
    nonseed_output="$output_dir/${base_name}${NON_SEED_SUFFIX}.tsv"

    # Step 1: Adding seed types to the input file
    echo "Adding seed types to $input_file..."
    python3 "$add_seeds_and_filtering_dir/add_seed_types.py" --ifile "$input_file" --ofile "$seed_types_output"
    if [ $? -ne 0 ]; then
        echo "Error in adding seed types. Check your script and input file."
        exit 1
    fi
    echo "Seed types added to $input_file. Output saved to $seed_types_output"

    # Step 2: Filtering interactions based on seed types
    echo "Filtering interactions for $seed_types_output..."
    python3 "$add_seeds_and_filtering_dir/filter_interactions.py" --ifile "$seed_types_output" --canonical_ofile "$canonical_output" --noncanonical_ofile "$noncanonical_output" --nonseed_ofile "$nonseed_output"
    if [ $? -ne 0 ]; then
        echo "Error in filtering interactions. Check your script and input file."
        exit 1
    fi
    echo "Interactions filtered for $seed_types_output. Canonical interactions saved to $canonical_output, noncanonical interactions saved to $noncanonical_output, and nonseed interactions saved to $nonseed_output"

done

echo "Adding seed_types and filtering interactions pipeline completed successfully."