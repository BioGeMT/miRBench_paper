#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=pp_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30

# parse command-line arguments
while getopts i:o: flag; do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

# check if required arguments are provided
if [ ! -d "$input_dir" ] || [ ! -d "$output_dir" ]; then
    echo "Usage: $0 -i input_dir -o output_dir"
    exit 1
fi

# define log file in the output directory
log_file="$output_dir/postprocess_1_genomic_region_annotation.log"

# redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# check for permission issues with log file
if [ $? -ne 0 ]; then
    echo "Error setting up log file. Check your permissions."
    exit 1
fi

# define paths to the directories where the scripts are located
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
filter_dir="$SCRIPT_DIR/../genomic_region_annotator_filtering"

# Ensembl release (fixed)
ENSEMBL_RELEASE=90

# define constants for suffixes with extensions
ANNOTATED_SUFFIX=".annotated"
SITE_SUMMARY_SUFFIX=".annotated_site_summary"
FILTERED_SUFFIX=".annotated.filtered"

# process each .tsv file in the input directory
for input_file in "$input_dir"/*.tsv; do
    if [ ! -f "$input_file" ]; then
        echo "No .tsv files found in $input_dir"
        exit 1
    fi

    base_name=$(basename "$input_file" .tsv)
    output_stem="$output_dir/${base_name}${ANNOTATED_SUFFIX}.tsv"
    summary_file="$output_dir/step2/${base_name}${SITE_SUMMARY_SUFFIX}.tsv"
    filtered_output="$output_dir/${base_name}${FILTERED_SUFFIX}.tsv"

    # Step 1: annotate
    echo "Running annotate for $input_file (Ensembl release $ENSEMBL_RELEASE) ..."
    genomic-region-annotator annotate \
        --input "$input_file" \
        --release "$ENSEMBL_RELEASE" \
        --output "$output_stem"
    if [ $? -ne 0 ]; then
        echo "Error in annotate step for $input_file"
        exit 1
    fi
    echo "Annotate completed for $input_file"

    # Locate Step 1 outputs (tool creates step1/ automatically)
    transcripts_file="$output_dir/step1/${base_name}${ANNOTATED_SUFFIX}_transcripts.tsv"
    matrix_file="$output_dir/step1/${base_name}${ANNOTATED_SUFFIX}_matrix.tsv"

    if [ ! -f "$transcripts_file" ] || [ ! -f "$matrix_file" ]; then
        echo "Error: Could not find Step 1 outputs for $base_name"
        echo "Expected: $output_dir/step1/${base_name}${ANNOTATED_SUFFIX}_transcripts.tsv and ${base_name}${ANNOTATED_SUFFIX}_matrix.tsv"
        exit 1
    fi

    # Step 2: summarize-sites
    echo "Running summarize-sites for $base_name ..."
    genomic-region-annotator summarize-sites \
        --transcripts "$transcripts_file" \
        --matrix "$matrix_file" \
        --policy clash_utr3_first \
        --dominance coverage \
        --report
    if [ $? -ne 0 ]; then
        echo "Error in summarize-sites step for $base_name"
        exit 1
    fi
    echo "Summarize-sites completed for $base_name"

    if [ ! -f "$summary_file" ]; then
        echo "Error: Could not find Step 2 site summary output for $base_name"
        echo "Expected: $summary_file"
        exit 1
    fi

    # Step 3: filter summarized output to downstream columns
    echo "Filtering summarized output for $summary_file ..."
    python3 "$filter_dir/genomic_region_annotator_filtering.py" --ifile "$summary_file" --ofile "$filtered_output"
    if [ $? -ne 0 ]; then
        echo "Error filtering summarized output for $summary_file"
        exit 1
    fi
    echo "Filtered output written to $filtered_output"
done

# Done
echo "postprocess_1 genomic region annotation pipeline completed successfully."
