
# Post-Processing Pipeline Script

This script is designed to process `.tsv` files through multiple stages including filtering, data splitting, and negative sample generation. 

## Requirements
- Python 3
- Run `conda env create --name <env_name> --file=post_process.yml`, then `conda activate <env_name>`
- Necessary Python scripts located in specified directories:
  - `filtering/filtering.py`
  - `make_neg_sets/make_neg_sets.py`

## Usage
```bash
./post_process.sh -i input_dir [-o output_dir] [-n intermediate_dir] [-t neg_ratios] [-r min_edit_distance]
```

## Options
- `-i` : Input directory containing `.tsv` files (required).
- `-o` : Output directory (optional, default: `output` in current directory).
- `-n` : Intermediate directory (optional, default: `intermediate` in current directory).
- `-t` : Comma-separated list of negative ratios (optional, default: `1,10,100`).
- `-r` : Minimum edit distance (optional, default: `3`).

## Script Directories
- `filtering_dir="../filtering"`
- `make_neg_sets_dir="../make_neg_sets"`

Ensure that the directories containing the scripts are correctly set relative to the location of the `pipeline.sh` script.

## Process Details
1. **Filtering**: The script filters the data using the `filtering.py` script.
   ```bash
   python3 filtering/filtering.py --ifile input_file --ofile filtered_file
   ```
   
2. **Deduplication**: The script deduplicates the filtered data using `awk`.
   ```bash
   awk -F'\t' 'NR==1{print $0} NR>1{if(!seen[$1$2]++){print}}' "$filtered_file" > "$deduplicated_file"
   ```

3. **Negative Sample Generation**: Negative samples are generated using the `make_neg_sets.py` script for each specified ratio.
   ```bash
   python3 make_neg_sets/make_neg_sets.py --ifile deduplicated_file --ofile neg_file --neg_ratio ratio --min_required_edit_distance min_required_edit_distance
   ```

4. **Data Splitting**: The data is split into train and test sets based on the `test` column using `awk`.
   ```bash
   awk -F'\t' 'NR==1{header=$0; print header > train_file; print header > test_file} NR>1{if($5=="False"){print > train_file} else {print > test_file}}' neg_file
   ```



## Log File
The script logs all output to a file named `pipeline.log` in the output directory.

## Example
```bash
./post_process.sh -i input_data -o output -n intermediate -t 1,10,100 -r 3
```

This example processes `.tsv` files in the `input_data` directory, outputs results to the `output` directory, and uses the `intermediate` directory for intermediate files. It generates negative samples with ratios 1, 10, and 100, and uses a minimum edit distance of 3.

## Notes
- Ensure the directories containing the Python scripts (`filtering`, `make_neg_sets`) are correctly set relative to the location of this script.
- The script creates the output and intermediate directories if they do not already exist.
