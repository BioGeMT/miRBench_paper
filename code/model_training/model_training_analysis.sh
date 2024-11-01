#!/bin/bash


# Parse command-line input
while getopts "i:" opt; do
  case $opt in
    i) input_file="$OPTARG" ;;
    *) echo "Usage: $0 -i <input_file>" >&2; exit 1 ;;
  esac
done

# Check if input_file is set, otherwise print usage and exit
if [ -z "$input_file" ]; then
  echo "Error: -i <input_file> is mandatory"
  echo "Usage: $0 -i <input_file>"
  exit 1
fi

# Make subsets of the dataset
echo "Creating dataset subsets..."
python make_subset.py --N 200 --dataset "$input_file" --output Manakov_subset_200.tsv
python make_subset.py --N 2000 --dataset "$input_file" --output Manakov_subset_2000.tsv
python make_subset.py --N 7720 --dataset "$input_file" --output Manakov_subset_7720.tsv
python make_subset.py --N 20000 --dataset "$input_file" --output Manakov_subset_20000.tsv
python make_subset.py --N 200000 --dataset "$input_file" --output Manakov_subset_200000.tsv

# Encode datasets
echo "Encoding datasets..."
python encode_dataset.py -i Manakov_subset_200.tsv -o Manakov_subset_200
python encode_dataset.py -i Manakov_subset_2000.tsv -o Manakov_subset_2000
python encode_dataset.py -i Manakov_subset_7720.tsv -o Manakov_subset_7720
python encode_dataset.py -i Manakov_subset_20000.tsv -o Manakov_subset_20000
python encode_dataset.py -i Manakov_subset_200000.tsv -o Manakov_subset_200000
python encode_dataset.py -i "$input_file" -o Manakov_full

# Train models
echo "Models training..."
python training.py --data Manakov_subset_200_dataset.npy --labels Manakov_subset_200_labels.npy --dataset_size 200 --ratio 1 --model CNN_Manakov_subset_200.keras
python training.py --data Manakov_subset_2000_dataset.npy --labels Manakov_subset_2000_labels.npy --dataset_size 2000 --ratio 1 --model CNN_Manakov_subset_2000.keras
python training.py --data Manakov_subset_7720_dataset.npy --labels Manakov_subset_7720_labels.npy --dataset_size 7720 --ratio 1 --model CNN_Manakov_subset_7720.keras
python training.py --data Manakov_subset_20000_dataset.npy --labels Manakov_subset_20000_labels.npy --dataset_size 20000 --ratio 1 --model CNN_Manakov_subset_20000.keras
python training.py --data Manakov_subset_200000_dataset.npy --labels Manakov_subset_200000_labels.npy --dataset_size 200000 --ratio 1 --model CNN_Manakov_subset_200000.keras
python training.py --data Manakov_full_dataset.npy --labels Manakov_full_labels.npy --dataset_size 2524246 --ratio 1 --model CNN_Manakov_full.keras
