#!/bin/bash

# Install gdown if not already installed
if ! command -v gdown &> /dev/null
then
    echo "gdown not found. Installing gdown..."
    pip install gdown
fi

# Define file IDs for the specific .tsv files you want to download
FILE_ID_1="15gQ4WShxUZjEjSwppc9cukfusanpbl6n"
FILE_ID_2="1fWD_5TVu1B6RWJ2aVm69M2TPbJ8Gh75Z"
FILE_ID_3="1_GYAxu5DzudQ5ka-zwqEa4friyMhF5vn"

# Create a directory to store the downloaded .tsv files
mkdir -p predictions

# Download each file by its ID
echo "Downloading AGO2_eCLIP_Manakov2022_100_CNN_predictions.tsv"
gdown "https://drive.google.com/uc?id=$FILE_ID_1" -O predictions/AGO2_eCLIP_Manakov2022_100_CNN_predictions.tsv

echo "Downloading AGO2_eCLIP_Manakov2022_10_CNN_predictions.tsv"
gdown "https://drive.google.com/uc?id=$FILE_ID_2" -O predictions/AGO2_eCLIP_Manakov2022_10_CNN_predictions.tsv

echo "Downloading AGO2_eCLIP_Manakov2022_1_CNN_predictions.tsv"
gdown "https://drive.google.com/uc?id=$FILE_ID_3" -O predictions/AGO2_eCLIP_Manakov2022_1_CNN_predictions.tsv

echo "Download completed! .tsv files are saved in the predictions/ directory."
