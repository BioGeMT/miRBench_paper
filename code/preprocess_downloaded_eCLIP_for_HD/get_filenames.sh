#!/bin/bash

## README
## To create a text file with a list of fastq filenames in a directory (excluding files ending with '_2.fastq.gz'), run:
## `bash get_filenames.sh <source directory> <output file>`

# Assign the first argument to the variable source_dir and the second argument to the variable output_file
source_dir="$1"
output_file="$2"

# Use the find command to list all files and redirect the output to the specified file
find "$source_dir" -type f -name "*.fastq.gz" -printf "%f\n" | grep -v "_2\.fastq\.gz$" > "$output_file"

echo "File list has been saved to $output_file"
