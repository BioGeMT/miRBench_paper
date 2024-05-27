# Preprocessing downloaded eCLIP data from Manakov 2022 for later processing in HybriDetector (HD)

## Preprocessing steps as per https://github.com/YeoLab/chim-eCLIP#extract-umi-and-trim (accessed by Stephanie Sammut on 21-May-2024);
1. Extract the 5' UMI from the reads to the read name
2. Trim the 3' adapters from the reads
3. Trim the 3' UMI from the reads

## To run the script, you will need:
1. To create a directory with all the files (or symlinks to the files) you wish to preprocess
    - Refer to create_symlinks.sh
2. To create a text file with a list of the names of the fastq files you wish to preprocess (one name per line, including ".fastq.gz" file extensions) 
    - Refer to get_filenames.sh
3. To choose a directory for the pipeline's outputs
4. To create and activate an environment with the following packages installed: 
    - umi_tools
    - cutadapt
    - gzip
5. To change the SBATCH directives in the script as required. Note that the script uses SLURM to run as a job array on a HPC, preprocessing multiple files in parallel. 
 
## From the command line, run:
`sbatch preprocess_eCLIP_for_HD.sh <data_directory> <filenames> <destination_directory>`
where 
`<data_directory>` is the directory containing the files (or symlinks to the files) to be preprocessed for HD 
`<filenames>` is the text file containing a list of file names to be preprocessed for HD 
`<destination_directory>` is a directory of choice to output the files to be run through HD 

## Note that the destination directory will have the following folder structure once all the pipeline has run:
dest_dir
    |- SAMPLE_NAME
    |   |- SAMPLE_NAME.pp.fq.gz
    |   |- logs/
    |- ANOTHER_SAMPLE_NAME
    
SAMPLE_NAME.pp.fq.gz is the file to be run through HybriDetector. `logs/` contains the reports from umi_tools and cutadapt, and the stdout of the script. The script creates a `temp/` folder containing temporary files, but it is deleted at the end of the script if the final output exists. 