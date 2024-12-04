#!/bin/bash


# Parse command-line input
while getopts i:o:n:s: opt; do
  case ${opt} in
    i) input_file=${OPTARG} ;;
    o) output_dir=${OPTARG};;
    n) intermediate_dir=${OPTARG};;
    s) IFS=',' read -r -a subsets <<< "${OPTARG}";;
    *) echo "Usage: $0 -i input_file [-o output_dir] [-n intermediate_dir] [-s subsets]
    Input is a tsv file with the dataset containing gene, noncodingRNA and label column. Subsets is a list of integers separated by commas.
    Output is a set of models trained on subsets of this dataset." >&2; exit 1 ;;
  esac
done

# Check if input_file is set, otherwise print usage and exit
if [ -z "$input_file" ]; then
  echo "Error: -i <input_file> is mandatory"
  echo "Usage: $0 -i input_file [-o output_dir] [-n intermediate_dir] [-s subsets]
    Input is a tsv file with the dataset containing gene, noncodingRNA and label column. Subsets is a list of integers separated by commas.
    Output is a set of models trained on subsets of this dataset."
  exit 1
fi

# set default values for subsets if not specified
default_subsets=(100 1000 3860 10000 100000 1262123)
subsets=( "${subsets[@]:-${default_subsets[@]}}" )

# define directories for output and intermediate files
output_dir="${output_dir:-$(pwd)/models}"
intermediate_dir="${intermediate_dir:-$(pwd)/intermediate}"

# create output and intermediate directories if they don't exist
mkdir -p "$intermediate_dir" "$output_dir"

# Make subsets, encode dataset and train models
for size in ${subsets[@]}; do
  echo "Creating subset of size $size..."

  python make_subset.py --N "$size" --dataset "$input_file" --output "${intermediate_dir}/subset_${size}.tsv"
  # Handle potential errors in the make_subset.py script - N is bigger than the dataset size
  if [ $? -ne 0 ]; then
        echo "Error making subset of size $size. Check if the dataset size is not smaller than $size."
        continue
  fi

  python encode_dataset.py -i "${intermediate_dir}/subset_${size}.tsv" -o "${intermediate_dir}/subset_${size}"

  echo "Training model"
  python training.py --data "${intermediate_dir}/subset_${size}_dataset.npy" --labels "${intermediate_dir}/subset_${size}_labels.npy" --dataset_size "$((size*2))" --ratio 1 --model "${output_dir}/CNN_subset_${size}.keras"
done