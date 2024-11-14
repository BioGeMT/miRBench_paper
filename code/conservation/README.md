# Addition of Conservation Scores

A Python script for processing genomic interval data and adding [conservation scores](./scores.md) from phyloP and phastCons BigWig files.

## Requirements

- Input file with 1-based genomic coordinates. The script adjusts them to 0-based to match the 0-based coordinate system in the following BigWig files. 
- phyloP and phastCons BigWig files need to be downloaded from:
  - https://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw 
  - https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.bw 
- Note that the conservation files used to create the Manakov datasets with conservation were downloaded on 15 September 2024 by Dimosthenis Tzimotoudis. 

- Required Python packages:
  - pandas
  - numpy
  - pyBigWig
  - argparse

## Usage

```bash
python add_conservation_scores.py --ifile input.tsv --ofile output.tsv --phyloP_path path/to/phyloP.bw --phastCons_path path/to/phastCons.bw
```

### Arguments

- `--ifile`: Path to input file with 1-based coordinates (TSV format). If not specified, reads from STDIN
- `--ofile`: Path to output file (TSV format). If not specified, writes to STDOUT
- `--phyloP_path`: Path to phyloP BigWig file with 0-based coordinates 
- `--phastCons_path`: Path to phastCons BigWig file with 0-based coordinates 

### Input File Format

The input file should be a tab-separated (TSV) file with the following required columns:
- chr (chromosome)
- start (1-based)
- end (1-based)
- strand
- gene (sequence)

### Output File Format

The script adds two new columns to the input data:
- `gene_phyloP`: Conservation scores from phyloP
- `gene_phastCons`: Conservation scores from phastCons

## Implementation Details

- Automatically adds 'chr' prefix if missing
- Converts "MT" to "chrM" for mitochondrial DNA
- Validates chromosome names against BigWig files
- Uses scores with 0-based input coordinates, converts input file genomic coordinates from 1- to 0- based for processing
- Checks interval boundaries and conservation score values
- Checks conservation scores length against gene sequence length
- Marks invalid entries with NaN values