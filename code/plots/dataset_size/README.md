# Dataset Size Plot 

The script visualizes the performance of machine a learning model across different dataset sizes using a logarithmic scale. It generates a plot showing the relationship between dataset size and model performance.

## Usage
Run the script from the command line with the following syntax:

```bash
python plot.py input.tsv output.png
```

### Arguments:
* `input.tsv`: Input TSV file containing model performance data
* `output.png`: Output PNG file for the generated plot

## Input File Format
The input TSV file should have the following structure:
* Header row with column names
* First column: Model/tool description
* Second column: Test set scores
* Last column: Dataset size

Example:
```
Tool / Dataset	Manakov test set	Manakov leftout	Hejret new test set	Klimentova new test set	Dataset Size
CNN_Manakov_subset_100	0.4841131007	0.4853891666	0.4866649963	0.4723581146	200
...
```

## Dependencies
The script requires the following Python libraries:
* numpy
* matplotlib
