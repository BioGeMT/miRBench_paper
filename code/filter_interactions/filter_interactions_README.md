# Filter interactions

## Overview

This script filters a file with annotated seed types into canonical, non-canonical, and non-seed interactions.

## Features

- Defines canonical, non-canonical, and non-seed interactions.
- Filters canonical, non-canonical, and non-seed interactions, and saves them into 3 distinct files.

## Requirements

- Python 3.8
- Required Python packages:
  - `pandas`
  - `argparse`

## Usage

Run the script with the following command:

```bash
python filter_interactions.py --ifile <input_file_with_seed_types> --canonical_ofile <output_file_with_canonical_interactions> --noncanonical_ofile <output_with_noncanonical_interactions> --nonseed <output_with_nonseed_interactions>
```

### Arguments

- `--ifile`: Path to input file with seed type annotations (columns: Seed6mer, Seed6merBulgeOrMismatch). This is output of add_seed_types.py script. 
- `--canonical_ofile`: Path to output file containing canonical (Seed6mer) interactions
- `--noncanonical_ofile`: Path to output file containing non-canonical (Seed6merBulgeOrMismatch but not Seed6mer) interactions. 
- `--nonseed_ofile`: Path to output file containing non-seed (No Seed6merBulgeOrMismatch) interactions.  

## Note

Note that in the miRBench package, Seed6merBulgeOrMismatch is the most loose seed type and therefore include all other seed types defined in the miRBench package, including the Seed6mer. 

