# Concatenate HybriDetector output file from different samples

This script concatenates all the *.unified_length_all_types_unique_high_confidence.tsv files in the input directory and saves the output to the output file. It is meant to concatenate some output files from the HybriDetector pipeline into a single file, in preparation for the post-processing pipeline. 

## Usage

`bash concat_HD_output.sh --input_dir <input_directory> --output_file <output_file>`

The script can also be submitted to a SLURM cluster using the following command:

`sbatch concat_HD_output.sh --input_dir <input_directory> --output_file <output_file>`

### Arguments

The script takes two arguments:

- `--input_dir`: The input directory containing the *.unified_length_all_types_unique_high_confidence.tsv files (or symlinks to such files)
- `--output_file`: The output file name to save the concatenated output

### Examples

`bash concat_HD_output.sh --input_dir input --output_file concatenated_HD_output.tsv`

`sbatch concat_HD_output.sh --input_dir input --output_file concatenated_HD_output.tsv`

## Notes

The script creates an output/ directory in the current directory. 
The <output_file> will be saved here. 
The script also logs all stdout and stderr to a log file which is also saved in the output/ directory. 

The script ensures that the header from the first file is included in the output file, while subsequent files have their headers removed to avoid duplication.

