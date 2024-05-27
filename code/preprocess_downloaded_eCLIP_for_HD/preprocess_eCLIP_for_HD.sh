#!/bin/bash
#SBATCH --account=ssamm10 
#SBATCH --job-name=preprocess_eCLIP_for_HD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --array=1-2%2 
#SBATCH --nice=500
#SBATCH --output=/dev/null

# Define paths
ECLIP_PATH="$1" 
ECLIP_FILE_LIST="$2" 
OUTPUT_PATH="$3" 

# Function to perform a step if output does not already exist
function run_step {
    local OUTPUT=$1
    local MSG_START=$2
    local MSG_END=$3
    local CMD=$4

    if [[ ! -f "$OUTPUT" ]]; then
        echo "$MSG_START" >> $LOGFILE
        eval $CMD && echo "$MSG_END" >> $LOGFILE
    else
        echo "$OUTPUT already exists, skipping." >> $LOGFILE
    fi
}

# Function to delete the temp folder if final out exists
function delete_temp {
    local OUTPUT=$1
    local MSG_START=$2
    local MSG_END=$3
    local CMD=$4

    if [[ -f "$OUTPUT" ]]; then
        echo "$MSG_START" >> $LOGFILE
        eval $CMD && echo "$MSG_END" >> $LOGFILE
    else
        echo "$OUTPUT does not exist, temp folder not deleted." >> $LOGFILE
    fi
}

# Read FASTQ filename from the list based on the SLURM array job index
ECLIP_FILE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $ECLIP_FILE_LIST)
ECLIP_FILE="$ECLIP_PATH/$ECLIP_FILE_NAME"
BASE_NAME=$(basename $ECLIP_FILE .fastq.gz)

# Create temporary output directory and its parents if it does not exist
if [[ ! -d $OUTPUT_PATH/$BASE_NAME/temp ]]; then
    mkdir -p $OUTPUT_PATH/$BASE_NAME/temp
fi

# Create logs directory and its parents if it does not exist
if [[ ! -d $OUTPUT_PATH/$BASE_NAME/logs ]]; then
   mkdir -p $OUTPUT_PATH/$BASE_NAME/logs
fi

# Define log file
LOGFILE="$OUTPUT_PATH/$BASE_NAME/logs/preprocess_eCLIP_for_HD.log"

# Step 1: Extract the 5' UMI from the reads to the read name
## Set the random seed to 1 to ensure reproducibility
## Set the barcode pattern to a 10-nt long UMI with each position being one of any of the four nucleotides

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

# Step 3: Trim the 3' UMI from the reads
## Set the number of bases to be removed from the 3' end of the reads to 10 (last 10 bases of the reads)
## Set the number of cores to use to 8

run_step "$OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.fq.gz" \
    "Extracting the 5' UMI from the reads to the read name for $BASE_NAME..." \
    "5' UMI extraction complete." \
    "umi_tools extract \
    --random-seed 1 \
    --stdin $ECLIP_FILE \
    --bc-pattern NNNNNNNNNN \
    --log $OUTPUT_PATH/$BASE_NAME/logs/$BASE_NAME.umi_tools.log \
    --stdout $OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.fq && gzip $OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.fq" && \
run_step "$OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.adapter.fq" \
    "Trimming the 3' adapters from the reads for $BASE_NAME..." \
    "Adapter trimming complete." \
    "cutadapt \
    -O 1 \
    -f fastq \
    --match-read-wildcards \
    --times 3 \
    -e 0.1 \
    --quality-cutoff 6 \
    -m 18 \
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
    $OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.fq.gz > $OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.adapter.fq 2> $OUTPUT_PATH/$BASE_NAME/logs/$BASE_NAME.3cutadapt.txt" && \
run_step "$OUTPUT_PATH/$BASE_NAME/$BASE_NAME.pp.fq.gz" \
    "Trimming the 3' UMI from the reads for $BASE_NAME..." \
    "3' UMI trimming complete." \
    "cutadapt \
    -u -10 \
    -j 8 \
    $OUTPUT_PATH/$BASE_NAME/temp/$BASE_NAME.umi.adapter.fq > $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.pp.fq 2> $OUTPUT_PATH/$BASE_NAME/logs/$BASE_NAME.5cutadapt.txt && gzip $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.pp.fq" && \
delete_temp "$OUTPUT_PATH/$BASE_NAME/$BASE_NAME.pp.fq.gz" \
    "Deleting temporary output files for $BASE_NAME..." \
    "Temporary output files deleted." \
    "rm -r $OUTPUT_PATH/$BASE_NAME/temp"

