
## Description

The `seed_utils.py` script contains utility functions for miRNA seed analysis. It provides functions to find seed matches between miRNA and target sequences, and to generate reverse complements of DNA sequences.

## Usage

This script is not meant to be run directly but to be imported and used by other scripts in the project: `violin.py`, `seeds.py`, and `families.py`.

To use these functions in another script, import them as follows:

```python
from seed_utils import reverse_complement, find_seed_match
```

## Dependencies

- Bio (Biopython)

Ensure that Biopython is installed in your Python environment.
