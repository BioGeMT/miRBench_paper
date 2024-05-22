#!/bin/bash
#SBATCH --account=ssamm10 
#SBATCH --job-name=preprocess_eCLIP_for_HD
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --array=1-47%10 
#SBATCH --nice=500
#SBATCH --output=preprocess_eCLIP_for_HD_%A_%a.out
#SBATCH --error=preprocess_eCLIP_for_HD_%A_%a.err 

# Define paths
ECLIP_PATH="$1" # /data/biogemt/miRNA_Binding_Site_Prediction/raw_data/eCLIP_Manakov_2022/data_preprocessed_for_HD/to_be_preprocessed-symlinks
ECLIP_FILE_LIST="$2" # to_be_preprocessed-filenames.txt (in current directory)
OUTPUT_PATH="$3" # /data/biogemt/miRNA_Binding_Site_Prediction/raw_data/eCLIP_Manakov_2022/data_preprocessed_for_HD

# Read FASTQ filename from the list based on the SLURM array job index
ECLIP_FILE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $ECLIP_FILE_LIST)
ECLIP_FILE="$ECLIP_PATH/$ECLIP_FILE_NAME"
BASE_NAME=$(basename $ECLIP_FILE .fastq.gz)

# Create output directory and temp output directory
mkdir -p $OUTPUT_PATH/$BASE_NAME
mkdir -p $OUTPUT_PATH/$BASE_NAME/temp

# Function to perform a step if output does not already exist
function run_step {
    local OUTPUT=$1
    local MSG_START=$2
    local MSG_END=$3
    local CMD=$4

    if [[ ! -f "$OUTPUT" ]]; then
        echo "$MSG_START"
        eval $CMD
        echo "$MSG_END"
    else
        echo "$OUTPUT already exists, skipping."
    fi
}

# Step 1: Extract the 5' UMI from the reads to the read name
## Set the random seed to 1 to ensure reproducibility
## Set the barcode pattern to a 10-nt long UMI with each position being one of any of the four nucleotides
echo "Extracting the 5' UMI from the reads to the read name for $BASE_NAME..."
run_step "$OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.fq.gz" \
    "Extracting the 5' UMI from the reads to the read name..." \
    "5' UMI extraction complete." \
    "umi_tools extract \
    --random-seed 1 \
    --stdin $ECLIP_FILE \
    --bc-pattern NNNNNNNNNN \
    --log $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.umi_tools.log \
    --stdout $OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.fq && gzip $OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.fq"

# Step 2: Trim the 3' adapters from the reads
## Set the minimum overlap length for adapter removal to 1
## Set the input format to fastq
## Allow IUPAC wildcards also in the reads
## Number of rounds of adapter matching per read set to 3
## Set the maximum error rate to 0.1
## Set the quality cutoff to 6 for trimming the 3' end of the reads
## Set the minimum length of the reads to 18
## Set the regular 3' adapters to be removed from the reads
## Set the number of cores to use to 8
echo "Trimming the 3' adapters from the reads for $BASE_NAME..."
run_step "$OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.adapter.fq" \
    "Trimming the 3' adapters from the reads..." \
    "Adapter trimming complete." \
    "cutadapt \
    -O 1 \ 
    -f fastq \
    --match-read-wildcards \
    --times 3 \
    -e 0.1 \
    --quality-cutoff 6 \
    -m 18 \
    -o $OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.adapter.fq \
    --json $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.adapters.cutadapt.json \
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
    $OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.fq.gz"

# Step 3: Trim the 3' UMI from the reads
## Set the number of bases to be removed from the 3' end of the reads to 10 (last 10 bases of the reads)
## Set the number of cores to use to 8
echo "Trimming the 3' UMI from the reads for $BASE_NAME..."
run_step "$OUTPUT_PATH/$BASE_NAME/$BASE_NAME.pp.fq.gz" \
    "Trimming the 3' UMI from the reads..." \
    "3' UMI trimming complete." \
    "cutadapt \
    -u -10 \
    -o $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.pp.fq \
    -j 8 \
    --json $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.umi.cutadapt.json \
    $OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.adapter.fq \
    && gzip $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.pp.fq" # pp for preprocessed

echo "Deleting temporary files for $BASE_NAME..."
rm -r $OUTPUT_PATH/$BASE_NAME/temp
    
