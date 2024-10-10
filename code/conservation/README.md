# Conservation Score Processing Scripts

This folder contains two Python scripts for adding and validating conservation scores:

1. `add_conservation_scores.py`: Adds conservation scores from BigWig files to data file with gene sequences.
2. `validate_conservation_scores.py`: Validates the conservation scores added by the first script.

## 1. add_conservation_scores.py

This script adds phyloP and phastCons conservation scores to data file with gene sequences. 

### Features:
- Reads gene data from a TSV file or standard input
- Processes data in blocks grouped by gene
- Adds phyloP and phastCons scores from BigWig files
- Handles special cases like mitochondrial chromosomes
- Outputs the data with added conservation scores

### Usage:
```
python add_conservation_scores.py --ifile <input_file> --ofile <output_file> --phyloP_path <phyloP_bigwig> --phastCons_path <phastCons_bigwig>
```

### Arguments:
- `--ifile`: Input TSV file (default: STDIN)
- `--ofile`: Output TSV file (default: STDOUT)
- `--phyloP_path`: Path to phyloP BigWig file
- `--phastCons_path`: Path to phastCons BigWig file

## 2. validate_conservation_scores.py

This script validates the conservation scores added by `add_conservation_scores.py`.

### Features:
- Reads data with conservation scores from a TSV file or standard input
- Validates that conservation scores are valid numbers
- Checks that the length of conservation scores matches the gene length
- Replaces invalid scores with NaN
- Outputs the data with validated conservation scores

### Usage:
```
python validate_conservation_scores.py --ifile <input_file> --ofile <output_file>
```

### Arguments:
- `--ifile`: Input TSV file (default: STDIN)
- `--ofile`: Output TSV file (default: STDOUT)

## Dependencies

Both scripts require the following Python libraries:
- pandas
- numpy
- pyBigWig (for `add_conservation_scores.py`)

## Notes

- Make sure your BigWig files (phyloP and phastCons) are compatible with your gene data file in terms of genome assembly and chromosome naming conventions.