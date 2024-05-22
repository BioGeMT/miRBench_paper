#!/bin/bash

# Assign the first argument to the variable source_dir and the second argument to the variable dest_dir
source_dir="$1"
dest_dir="$2"

# Create a symbolic link to the files in the source directory and its subfolders in the destination directory
find "$source_dir" -type f -name "*.fastq.gz" ! -name "*_2.fastq.gz" -exec ln -s {} "$dest_dir" \;

