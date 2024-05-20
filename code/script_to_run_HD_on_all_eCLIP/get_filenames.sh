#!/bin/bash

# Assign the first argument to the variable directory and the second argument to the variable output_file
directory="$1"
output_file="$2"

# Use the find command to list all files and redirect the output to the specified file
find "$directory" -type f -name "*.fastq.gz" | grep -v "_2\.fastq\.gz$" > "$output_file"

echo "File list has been saved to $output_file"
