
# Post-Processing Pipeline Script

This script is designed to process `.tsv` files through multiple stages including filtering, family assignment, data splitting, and negative sample generation. 

## Requirements
- Python 3
- Necessary Python scripts located in specified directories:
  - `filtering/filtering.py`
  - `family_assign/family_assign.py`
  - `make_neg_sets/make_neg_sets.py`
- `wget` for downloading the `mature.fa` file

## Usage
\`\`\`bash
./pipeline.sh -i input_dir [-o output_dir] [-n intermediate_dir] [-t neg_ratios] [-r min_edit_distance]
\`\`\`

## Options
- \`-i\` : Input directory containing `.tsv` files (required).
- \`-o\` : Output directory (optional, default: \`output_data\` in current directory).
- \`-n\` : Intermediate directory (optional, default: \`intermediate_data\` in current directory).
- \`-t\` : Comma-separated list of negative ratios (optional, default: \`1,10,100\`).
- \`-r\` : Minimum edit distance (optional, default: \`3\`).

## Script Directories
- \`filtering_dir="../filtering"\`
- \`family_assign_dir="../family_assign"\`
- \`make_neg_sets_dir="../make_neg_sets"\`
- \`pre_process_dir="../pre_process"\`

Ensure that the directories containing the scripts are correctly set relative to the location of the \`pipeline.sh\` script.

## Process Details
1. **Filtering**: The script filters the data using the \`filtering.py\` script.
   \`\`\`bash
   python3 filtering/filtering.py --ifile input_file --ofile filtered_file
   \`\`\`
   
2. **Family Assignment**: The script assigns families using the \`family_assign.py\` script and the \`mature.fa\` file.
   \`\`\`bash
   python3 family_assign/family_assign.py --ifile filtered_file --mature mature_file --ofile family_assigned_file
   \`\`\`

3. **Data Splitting**: The data is split into train and test sets based on the \`test\` column using \`awk\`.
   \`\`\`bash
   awk -F'\t' 'NR==1{header=$0; print header > train_file; print header > test_file} NR>1{if($5=="False"){print > train_file} else {print > test_file}}' family_assigned_file
   \`\`\`

4. **Negative Sample Generation**: Negative samples are generated for both train and test sets using the \`make_neg_sets.py\` script for each specified ratio.
   \`\`\`bash
   python3 make_neg_sets/make_neg_sets.py --ifile train_file --ofile output_file --neg_ratio ratio --min_edit_distance min_edit_distance
   \`\`\`

## Log File
The script logs all output to a file named \`pipeline.log\` in the output directory.

## Example
\`\`\`bash
./pipeline.sh -i data/input -o data/output -n data/intermediate -t 1,10,100 -r 3
\`\`\`

This example processes `.tsv` files in the \`data/input\` directory, outputs results to the \`data/output\` directory, and uses the \`data/intermediate\` directory for intermediate files. It generates negative samples with ratios 1, 10, and 100, and uses a minimum edit distance of 3.

## Notes
- Ensure the directories containing the Python scripts (\`filtering\`, \`family_assign\`, \`make_neg_sets\`, \`pre_process\`) are correctly set relative to the location of this script.
- The \`mature.fa\` file will be downloaded if not already present in the \`mature_data\` directory.
- The script creates the output and intermediate directories if they do not already exist.
