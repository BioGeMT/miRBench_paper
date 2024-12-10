#!/bin/bash

# Install gdown if not already installed
if ! command -v gdown &> /dev/null
then
    echo "gdown not found. Installing gdown..."
    pip install gdown
fi

download_dir="predictions"
plot_output_dir="output"
FILE_ID_1="15gQ4WShxUZjEjSwppc9cukfusanpbl6n"

mkdir -p $download_dir
mkdir -p $plot_output_dir

echo "Downloading AGO2_eCLIP_Manakov2022_100_CNN_predictions.tsv"
gdown "https://drive.google.com/uc?id=$FILE_ID_1" -O $download_dir/AGO2_eCLIP_Manakov2022_100_CNN_predictions.tsv

echo "Download completed! .tsv files are saved in the $download_dir directory."
