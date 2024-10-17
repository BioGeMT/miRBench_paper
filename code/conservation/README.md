# Addition of Conservation Scores

A Python script for processing genomic interval data and adding conservation scores from phyloP and phastCons BigWig files.

## Requirements

- Required Python packages:
  - pandas
  - numpy
  - pyBigWig
  - argparse

## Usage

```bash
python script.py --ifile input.tsv --ofile output.tsv --phyloP_path path/to/phyloP.bw --phastCons_path path/to/phastCons.bw
```

### Arguments

- `--ifile`: Input file path (TSV format). If not specified, reads from STDIN
- `--ofile`: Output file path (TSV format). If not specified, writes to STDOUT
- `--phyloP_path`: Path to the phyloP BigWig file
- `--phastCons_path`: Path to the phastCons BigWig file

### Input File Format

The input file should be a tab-separated (TSV) file with the following required columns:
- chr (chromosome)
- start
- end
- strand
- gene (sequence)

### Output

The script adds two new columns to the input data:
- `gene_phyloP`: Conservation scores from phyloP
- `gene_phastCons`: Conservation scores from phastCons

## Implementation Details

- Automatically adds 'chr' prefix if missing
- Converts "MT" to "chrM" for mitochondrial DNA
- Validates chromosome names against BigWig files
- Uses 1-based input coordinates, converts to 0-based for processing
- Checks interval boundaries and conservation score values
- Checks conservation scores length against gene sequence length
- Marks invalid entries with NaN values

## Example

```bash
python add_conservation_scores.py \
  --ifile genes.tsv \
  --ofile genes_with_conservation.tsv \
  --phyloP_path hg38.phyloP100way.bw \
  --phastCons_path hg38.phastCons100way.bw
```