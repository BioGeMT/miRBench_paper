# Preprocessing downloaded eCLIP data

Preprocessing steps done as per https://github.com/YeoLab/chim-eCLIP#extract-umi-and-trim (accessed by Stephanie Sammut on 21-May-2024);
1. Extract the 5' UMI from the reads to the read name
2. Trim the 3' adapters from the reads
3. Trim the 3' UMI from the reads

Make sure you have the following installed: 
- umi_tools
- cutadapt
- gzip

Change the SBATCH directives in preprocess_eCLIP_for_HD.sh as required. 

From the command line, run:

`sbatch preprocess_eCLIP_for_HD.sh <data_directory> <destination_directory>`

or, if not using the SLURM job scheduler, run:

`bash preprocess_eCLIP_for_HD.sh <data_directory> <destination_directory>`