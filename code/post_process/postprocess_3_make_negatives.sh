#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=pp_3
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
log_file="$output_dir/postprocess_3_make_negatives.log"

# redirect all output to the log file
exec > >(tee -a "$log_file") 2>&1

# check for permission issues with log file
if [ $? -ne 0 ]; then
    echo "Error setting up log file. Check your permissions."
    exit 1
fi

# define paths to the directories where the scripts are located
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
clustering_dir="$SCRIPT_DIR/../clustering"
make_negs_dir="$SCRIPT_DIR/../make_neg_sets"

# define constants for suffixes with extensions
GENE_ID_LOOKUP_SUFFIX=".gene_id_lookup"
CLUSTERING_OUTPUT_SUFFIX=".gene_clusters"
CLUSTERS_ADDED_SUFFIX=".gene_clusters_added"
NEG_SUFFIX=".negatives"

for input_file in "$input_dir"/*.tsv; do
    # process the input file
    base_name=$(basename "$input_file" .tsv)
    fasta_file="$intermediate_dir/${base_name}.fasta"
    gene_lookup_file="$intermediate_dir/${base_name}${GENE_ID_LOOKUP_SUFFIX}.tsv"
    clustering_output="$intermediate_dir/${base_name}${CLUSTERING_OUTPUT_SUFFIX}.csv"
    input_file_with_clusters="$intermediate_dir/${base_name}${CLUSTERS_ADDED_SUFFIX}.tsv"
    neg_output="$output_dir/${base_name}${CLUSTERS_ADDED_SUFFIX}${NEG_SUFFIX}.tsv"

    # Step 1: Generating FASTA file
    echo "Generating FASTA file for $input_file..."
    python3 "$clustering_dir/gene_fasta.py" --input "$input_file" --output "$fasta_file" --lookup "$gene_lookup_file"
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
    python3 "$clustering_dir/map_gene_clusters.py" --cluster_csv "$clustering_output" --dataset_tsv "$input_file" --lookup_tsv "$gene_lookup_file" --output_tsv "$input_file_with_clusters"
    if [ $? -ne 0 ]; then
        echo "Error in mapping clusters to input file. Check your script and input file."
        exit 1
    fi
    echo "Clusters mapped to $input_file. Output saved to $input_file_with_clusters"

    # Step 4: Make negatives directly from the clustered TSV
    echo "Generating negatives for $input_file_with_clusters..."
    python3 "$make_negs_dir/make_neg_sets.py" --ifile "$input_file_with_clusters" --ofile "$neg_output"
    if [ $? -ne 0 ]; then
        echo "Error in generating negative samples. Check your script and input file."
        exit 1
    fi
    echo "File with negative samples saved to $neg_output"
done

echo "postprocess_3 make negatives pipeline completed successfully."
