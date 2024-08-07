# Benchmark all tools on all datasets

## Instructions

### Create and activate a new python environment with python 3.8 

`python3.8 -m venv mirbench_env`

`source mirbench_env/bin/activate`

Note: miRBench in a conda environment still needs be tested - ran into issues

### Install miRBench

`pip install git+https://github.com/katarinagresova/miRBench.git`

### Run benchmark_all.py

Submit run_benchmark_all.sh to an HPC cluster via sbatch

`sbatch run_benchmark_all.sh -o <output_directory> -d <download_directory` 


