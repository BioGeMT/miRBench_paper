# Concatenate HD Post-process Pipeline Outputs

This script concatenates the output files from the HD post-processing pipeline into a single file per type and size. 

## Usage

`bash concat_post_process_output.sh -i <input_directory> -o <output_directory>`

The script can also be submitted to a SLURM cluster using the following command:

`sbatch concat_post_process_output.sh -i <input_directory> -o <output_directory>`

## Examples

`bash concat_post_process_output.sh -i input -o output`

`sbatch concat_post_process_output.sh -i input -o output`

### Arguments

- `-i`: Input directory containing all output files from the HD post processing pipeline
- `-o`: Output directory where the concatenated output files will be saved
  
## Description

- The script loops through each dataset type (train, test) and size (1, 10, 100), and concatenates multiple files with the same type and size into a single output file.
- The header line is taken from the first file, and the content is appended from all matching files, skipping the header line.
- Sample input and output files are available in the input/ and output/ dirs, respectively. 
- For the given types and sizes, 6 output files should be produced. 
- The output files are named based on the dataset type and size, e.g. AGO2_eCLIP_Manakov2022_1_train_dataset.tsv
- The script also logs all output to a log file in the output directory.

