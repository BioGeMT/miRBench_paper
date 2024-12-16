#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=pp_2
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
log_file="$output_dir/postprocess_2_make_negatives.log"

# redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# check for permission issues with log file
if [ $? -ne 0 ]; then
    echo "Error setting up log file. Check your permissions."
    exit 1
fi

# define paths to the directories where the scripts are located
clustering_dir="../clustering"
make_negs_dir="../make_neg_sets"

# define constants for suffixes with extensions
CLUSTERING_OUTPUT_SUFFIX=".gene_clusters"
CLUSTERS_ADDED_SUFFIX=".gene_clusters_added"
SORTED_SUFFIX=".mirfam_sorted"
NEG_SUFFIX=".negatives"

for input_file in "$input_dir"/*.tsv; do
    # process the input file
    base_name=$(basename "$input_file" .tsv)
    fasta_file="$intermediate_dir/${base_name}.fasta"
    clustering_output="$intermediate_dir/${base_name}${CLUSTERING_OUTPUT_SUFFIX}.csv"
    input_file_with_clusters="$intermediate_dir/${base_name}${CLUSTERS_ADDED_SUFFIX}.tsv"
    mirfam_sorted_file="$intermediate_dir/${base_name}${CLUSTERS_ADDED_SUFFIX}${SORTED_SUFFIX}.tsv"
    neg_output="$output_dir/${base_name}${CLUSTERS_ADDED_SUFFIX}${SORTED_SUFFIX}${NEG_SUFFIX}.tsv"

    # Step 1: Generating FASTA file
    echo "Generating FASTA file for $input_file..."
    python3 "$clustering_dir/gene_fasta.py" --input "$input_file" --output "$fasta_file"
    if [ $? -ne 0 ]; then
        echo "Error in generating FASTA file for $input_file. Check your script and input file."
        exit 1
    fi
    echo "FASTA file generated for $input_file. Output saved to $fasta_file"

    # Step 2: Performing sequence clustering on gene sequences in the generated FASTA file
    echo "Running gene sequence clustering for $fasta_file..."
    Rscript "$clustering_dir/clustering.R" "$fasta_file" "$clustering_output"
    if [ $? -ne 0 ]; then
        echo "Error in running gene sequence clustering. Check your script and input file."
        exit 1
    fi
    echo "Gene sequence clustering completed. Output saved to $clustering_output"

    # Step 3: Mapping clusters to input file
    echo "Mapping clusters to $input_file..."
    python3 "$clustering_dir/map_gene_clusters.py" --cluster_csv "$clustering_output" --dataset_tsv "$input_file" --output_tsv "$input_file_with_clusters"
    if [ $? -ne 0 ]; then
        echo "Error in mapping clusters to input file. Check your script and input file."
        exit 1
    fi
    echo "Clusters mapped to $input_file. Output saved to $input_file_with_clusters"

    # Step 4: Sort the file based on the noncodingRNA_fam column in preparation for negative sample generation
    echo "Sorting the input file with added clusters based on the noncodingRNA_fam column..."

    # Find the column number of the "noncodingRNA_fam" column
    column_number=$(head -n 1 "$input_file_with_clusters" | tr '\t' '\n' | nl -v 0 | grep "noncodingRNA_fam" | awk '{print $1}')

    # If the column number is found, sort the file by that column
    if [ -n "$column_number" ]; then
        (head -n 1 "$input_file_with_clusters" && tail -n +2 "$input_file_with_clusters" | sort -k "$column_number") > "${mirfam_sorted_file}"
        echo "Input file with added clusters sorted by the 'noncodingRNA_fam' column. Output saved to $mirfam_sorted_file"
    else
        echo "Error: 'noncodingRNA_fam' column not found in $mirfam_sorted_file"
        exit 1
    fi

    # Step 5: Make negatives
    echo "Generating negatives for $mirfam_sorted_file..."
    python3 "$make_negs_dir/make_neg_sets.py" --ifile "$mirfam_sorted_file" --ofile "$neg_output"
    if [ $? -ne 0 ]; then
        echo "Error in generating negative samples. Check your script and input file."
        exit 1
    fi
    echo "File with negative samples saved to $neg_output"

done

echo "Negative samples generation successfully completed for $input_file"