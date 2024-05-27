#!/bin/bash

## README
## To create symlinks in a destination directory, to files in a source directory (excluding files ending with '_2.fastq.gz'), run:
## `bash create_symlinks.sh <source directory> <destination directory>`

# Assign the first argument to the variable source_dir and the second argument to the variable dest_dir
source_dir="$1"
dest_dir="$2"

# Create a symbolic link to the files in the source directory and its subfolders in the destination directory
find "$source_dir" -type f -name "*.fastq.gz" ! -name "*_2.fastq.gz" -exec ln -s {} "$dest_dir" \;

