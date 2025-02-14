# Seed Type Annotation

## Overview

This script annotates seed types using the miRBench package.

## Features

- Indicates the presence of Seed6mer and Seed6merBulgeOrMismatch per noncodingRNA:gene pair, using the **miRBench** seed encoders and predictors.

## Requirements

- Python 3.8
- Required Python packages:
  - `miRBench`
  - `pandas`
  - `argparse`

## Usage

Run the script with the following command:

```bash
python add_seed_types.py --ifile <input_file> --ofile <output_file>
```

### Arguments

- `--ifile` (required): Path to input file with `noncodingRNA` and `gene` columns
- `--ofile` (required): Path to output file with added seed types (columns: `Seed6mer`, `Seed6merBulgeOrMismatch`)


