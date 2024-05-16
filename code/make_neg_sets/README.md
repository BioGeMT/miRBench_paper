input file sample

seq.g	seq.m	noncodingRNA_fam	feature	test	label
TCCAAATTAGAAGGCCGGCCCCGTCCACATTTTGCACAGTGCCTTTACAG	TGTGCAAATCCATGCAAAACTGA	mir-19	three_prime_utr	True	1
TATTTACTATAATGTTAGCTTACAAGCTGGGAATATAAGTGCATTAACGG	TATTGCACTTGTCCCGGCCTGT	mir-25	three_prime_utr	False	1
ATATTCTCAGGCCGCAAGTGCAATGCCTGAGGGGATCAGGCTTTTCTACT	TATTGCACTTGTCCCGGCCTGT	mir-25	three_prime_utr	False	1
ACTAAGGAACTGCAGCATTTGCACAGGGGAGGGGGGTGCCCTCCTTCCTA	TGTGCAAATCCATGCAAAACTGA	mir-19	three_prime_utr	False	1

--

I have this type of file.
These are my 'positive label (1)' samples.

I want to create N (e.g. 100) negative label (0) samples for each of the positives. 
I will do that by associating the same seq.g with a different seq.m randomly selected in this file. 
However, I want that seq.m to be a certain 'edit distance' away from the original one. Levenshtein edit distance.

How to run:

python make_neg_sets.py \
    --ifile pairs.tab \
    --ofile pairs.posneg.tab \
    --neg_ratio 100 \
    --min_edit_distance 3 \

## Prerequisites

Make sure you have the following Python packages installed:

- `pandas`
- `python-Levenshtein`

You can install these packages using pip:

```sh
pip install pandas python-Levenshtein