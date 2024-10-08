
# Post-Processing Pipeline Script

This script is designed to process `.tsv` files through multiple stages including filtering, family assignment, negative sample generation (at different class imbalance ratios), data splitting into test and train sets, and addition of conservation scores. 

## Requirements
- Python 3
- Run `conda env create --name <env_name> --file=post_process.yml`, then `conda activate <env_name>`
- Necessary Python scripts located in specified directories relative to the post_process.sh script:
  - `../filtering/filtering.py`
  - `../family_assign/family_assign.py`
  - `../make_neg_sets/make_neg_sets.py`
  - `../conservation/add_conservation_scores.py`
  - `../conservation/validate_conservation_scores.py`
- `wget` for downloading the `mature.fa` file

## Usage
```bash
./post_process.sh -i input_dir -p phyloP_path -c phastCons_path [-o output_dir] [-n intermediate_dir] [-t neg_ratios] [-r min_edit_distance]
```

## Options
- `-i` : Input directory containing `.tsv` files (required).
- `-p` : Path to phyloP BigWig file (required). 
- `-c` : Path to phastCons BigWig file (required). 
- `-o` : Output directory (optional, default: `output` in current directory).
- `-n` : Intermediate directory (optional, default: `intermediate` in current directory).
- `-t` : Comma-separated list of negative ratios (optional, default: `1,10,100`).
- `-r` : Minimum edit distance (optional, default: `3`).

## Process Details
1. **Filtering**: The script filters the data for the noncodingRNA 'mirna' and for specific columns.
   ```bash
   python3 ../filtering/filtering.py --ifile input_file --ofile filtered_file
   ```
   
2. **Deduplication**: The command deduplicates the filtered data based on concatenation of the gene and noncodingRNA sequences using `awk`.
   ```bash
   awk -F'\t' 'NR==1{print $0} NR>1{if(!seen[$1$2]++){print}}' "$filtered_file" > "$deduplicated_file"
   ```

3. **Family Assignment**: The script downloads the `mature.fa` file and assigns a family from it to noncodingRNA that are not assigned a noncodingRNA_fam. 
   ```bash
   python3 ../family_assign/family_assign.py --ifile filtered_file --mature mature_file --ofile family_assigned_file
   ```

4. **Sorting by Gene Sequence**: The command sorts the gene column alphabetically in preparation for the next step - negative sample generation. 
```bash
(head -n 1 "$family_assigned_file" && tail -n +2 "$family_assigned_file" | sort -k 1) > "${family_assigned_file_sorted}"
```

5. **Negative Sample Generation**: The script generates negative examples for each specified ratio.
   ```bash
   python3 ../make_neg_sets/make_neg_sets.py --ifile family_assigned_file --ofile neg_file --neg_ratio ratio --min_required_edit_distance min_required_edit_distance
   ```

6. **Data Splitting**: The data is split into train and test sets based on the `test` column using `awk`.
   ```bash
   awk -F'\t' 'NR==1{header=$0; print header > "'"$train_file"'"; print header > "'"$test_file"'"} NR>1{if($5=="False"){print > "'"$train_file"'"} else {print > "'"$test_file"'"}}' "$neg_file"
   ```

7. **Dropping the fifth (test) column**: This command removes the 'test' column used for the data splits as it is no longer needed. 
   ```bash
   awk -F'\t' 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) if(i!=5) printf "%s%s", $i, (i==NF?"\n":OFS)}' "$file" > "${file}_tmp" && mv "${file}_tmp" "$file"
   ```

8. **Addition of Conservation scores**: This script adds phyloP and phastCons conservation scores to the data. 
   ```bash
   python3 "$conservation_dir/add_conservation_scores.py" --ifile "$file" --phyloP "$phyloP_path" --phastCons "$phastCons_path" --ofile "$conservation_file"
   ```

9. **Validation of Conservation scores**: This script validates conservaton scores values and length, relative to the gene sequence length. 
   ```bash
   python3 "$conservation_dir/validate_conservation_scores.py" --ifile "$conservation_file" --ofile "$cleaned_conservation_file"
   ```

## Log File
The script logs all output to a file named `pipeline.log` in the output directory.

## Notes
- The `mature.fa` file will only be downloaded if not already present in the `mature` directory.
