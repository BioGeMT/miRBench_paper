#!/bin/bash
#SBATCH --account=ssamm10
#SBATCH --job-name=preprocess_eCLIP_for_HD
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --array=1-47%6 
#SBATCH --nice=500
#SBATCH --output=preprocess_eCLIP_for_HD_%A_%a.out
#SBATCH --error=preprocess_eCLIP_for_HD_%A_%a.err 

# Assign the first argument to the variable data_dir, the second argument to the variable filename_file, and the third to the variable dest_dir
data_dir="$1" # /data/biogemt/miRNA_Binding_Site_Prediction/raw_data/eCLIP_Manakov_2022/data_symlinks
filename_file="$2" # /data/biogemt/miRNA_Binding_Site_Prediction/raw_data/eCLIP_Manakov_2022/filenames.txt
dest_dir="$3"  # /data/biogemt/miRNA_Binding_Site_Prediction/raw_data/eCLIP_Manakov_2022/data_preprocessed_for_HD

# Define paths
ECLIP_PATH=$data_dir
ECLIP_FILE_LIST=$filename_file
OUTPUT_PATH=$dest_dir

# Read FASTQ filename from the list based on the SLURM array job index
ECLIP_FILE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $ECLIP_FILE_LIST)
ECLIP_FILE="$ECLIP_PATH/$ECLIP_FILE_NAME"
BASE_NAME=$(basename $ECLIP_FILE .fastq.gz)

# Create temp output directory and move into it
mkdir -p $OUTPUT_PATH/temp_out
cd $OUTPUT_PATH/temp_out

# For each file in the data directory (as per https://github.com/YeoLab/chim-eCLIP?tab=readme-ov-file#extract-umi-and-trim; accessed by Stephanie Sammut on 21-May-2024):
# 1. Extract the 5' UMI from the reads to the read name
# 2. Trim the 3' adapters from the reads
# 3. Trim the 3' UMI from the reads
for file in "$data_dir"/*.fastq.gz
do
    # Get the base name of the file
    base_name=$(basename "$file" .fastq.gz)
    
    # Extract the 5' UMI from the reads to the read name
    ## Set the random seed to 1 to ensure reproducibility
    ## Read the input file
    ## Set the barcode pattern to a 10-nt long UMI with each position being one of any of the four nucleotides
    ## Write the output to a file
    umi_tools extract \ 
    --random-seed 1 \ 
    --stdin "$file" \ 
    --bc-pattern NNNNNNNNNN \ 
    --log ${file}.log \ # Log the output
    --stdout ${base_name}.umi.fq 

    gzip ${base_name}.umi.fq

    # Trim the 3' adapters from the reads
    ## Set the minimum overlap length for adapter removal to 1
    ## Set the input format to fastq
    ## Allow IUPAC wildcards also in the reads
    ## Number of rounds of adapter matching per read set to 3
    ## Set the maximum error rate to 0.1
    ## Set the quality cutoff to 6 for trimming the 3' end of the reads
    ## Set the minimum length of the reads to 18
    ## Write the output to a file
    ## Set the regular 3' adapters to be removed from the reads
    ## Set the number of cores to use to 8
    ## Read the input file
    cutadapt \
    -O 1 \ 
    -f fastq \
    --match-read-wildcards \
    --times 3 \
    -e 0.1 \
    --quality-cutoff 6 \
    -m 18 \
    -o ${base_name}.umi.fqTr.fq \
    -a AGATCGGAAG \
    -a GATCGGAAGA \
    -a ATCGGAAGAG \
    -a TCGGAAGAGC \
    -a CGGAAGAGCA \
    -a GGAAGAGCAC \
    -a GAAGAGCACA \
    -a AAGAGCACAC \
    -a AGAGCACACG \
    -a GAGCACACGT \
    -a AGCACACGTC \
    -a GCACACGTCT \
    -a CACACGTCTG \
    -a ACACGTCTGA \
    -a CACGTCTGAA \
    -a ACGTCTGAAC \
    -a CGTCTGAACT \
    -a GTCTGAACTC \
    -a TCTGAACTCC \
    -a CTGAACTCCA \
    -a TGAACTCCAG \
    -a GAACTCCAGT \
    -a AACTCCAGTC \
    -a ACTCCAGTCA \
    -j 8 \
    ${base_name}.umi.fq.gz

    # Trim the 3' UMI from the reads
    # Set the number of bases to be removed from the 3' end of the reads to 10 (last 10 bases of the reads)
    # Write the output to a file
    # Set the number of cores to use to 8
    # Read the input file
    cutadapt \
    -u -10 \
    -o ${base_name}.umi.fqTrTr.fq.fastq \
    -j 8 \
    ${base_name}.umi.fqTr.fq
    
    gzip ${base_name}.umi.fqTrTr.fq.fastq

    # Move the preprocessed file to the destination directory
    mv ${base_name}.umi.fqTrTr.fq.fastq.gz $dest_dir/$base_name.umi.trm.fastq.gz
    
    # Move back to the destination directory and remove the intermediate files directory
    cd ..
    rm -r intermediate_files
    
    # Indicate that the preprocessing of the file is done
    echo "Preprocessing of $file is done."
done




