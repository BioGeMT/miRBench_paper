# Concatenate HybriDetector output file from different samples

This script concatenates all the *.unified_length_all_types_unique_high_confidence.tsv files in the input directory and saves the output to the output directory. It is meant to concatenate some output files from the HybriDetector pipeline into a single file, in preparation for the post-processing pipeline. 

## Usage

`bash concat_HD_output.sh --input_dir <input_directory> --output_dir <output_dir> --output_prefix <output_prefix>`

The script can also be submitted to a SLURM cluster using the following command:

`sbatch concat_HD_output.sh --input_dir <input_directory> --output_dir <output_dir> --output_prefix <output_prefix>`

### Arguments

The script takes two arguments:

- `--input_dir`: The input directory containing the *.unified_length_all_types_unique_high_confidence.tsv files (or symlinks to such files)
- `--output_dir`: The output directory containing the output file and log file
- `--output_prefix`: The prefix to the output filename which will be followed by .unified_length_all_types_unique_high_confidence.tsv (default: concat)

### Examples

`bash concat_HD_output.sh --input_dir input --output_dir output --output_prefix concat`

`sbatch concat_HD_output.sh --input_dir input --output_dir output --output_prefix concat`

## Notes

The script creates the output directory in the current directory if it does not already exist. 
The output file will be saved here as `<prefix>.unified_length_all_types_unique_high_confidence.tsv`
The script also logs all stdout and stderr to a log file which is also saved in the output directory. 

The script ensures that the header from the first file is included in the output file, while subsequent files have their headers removed to avoid duplication.

