# Preprocessing downloaded eCLIP data

## To create a text file with a list of fastq filenames in a directory (excluding files ending with '_2.fastq.gz'), run:
`bash get_filenames.sh <source directory> <output file>`

## To create symlinks in a destination directory, to files in a source directory (excluding files ending with '_2.fastq.gz'), run:
`bash create_symlinks.sh <source directory> <destination directory>`

## Preprocessing steps done as per https://github.com/YeoLab/chim-eCLIP#extract-umi-and-trim (accessed by Stephanie Sammut on 21-May-2024);
1. Extract the 5' UMI from the reads to the read name
2. Trim the 3' adapters from the reads
3. Trim the 3' UMI from the reads

Make sure you have the following installed: 
- umi_tools
- cutadapt
- gzip

The script uses SLURM to run as a job array, preprocessing multiple files in parallel. 

Change the SBATCH directives in preprocess_eCLIP_for_HD.sh as required. 

From the command line, run:
`sbatch preprocess_eCLIP_for_HD.sh <data_directory> <filenames> <destination_directory>`

where 
`<data_directory>` is the directory containing symlinks to the files to be preprocessed for HD (produced by create_symlinks.sh)
`<filenames>` is the text file containing a list of file names to be preprocessed for HD (produced by get_filenames.sh)
`<destination_directory>` is a directory of choice to output the files to be run through HD 

e.g. 

sbatch preprocess_eCLIP_for_HD.sh /data/biogemt/miRNA_Binding_Site_Prediction/raw_data/eCLIP_Manakov_2022/data_preprocessed_for_HD/to_be_preprocessed-symlinks to_be_preprocessed-filenames.txt /data/biogemt/miRNA_Binding_Site_Prediction/raw_data/eCLIP_Manakov_2022/data_preprocessed_for_HD

Note that the destination directory will have the following folder structure:
dest_dir
    |- SAMPLE_NAME
    |   |- SAMPLE_NAME.umi_tools.log
    |   |- SAMPLE_NAME.pp.fq.gz
    |   |- temp/
    |- ANOTHER_SAMPLE_NAME

SAMPLE_NAME.pp.fq.gz is the file to be run through HybriDetector. The rest are log files from the commands run on that sample. `temp/` has all the intermediate files and can be deleted. 