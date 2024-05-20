#!/bin/bash
#SBATCH --account=ssamm10
#SBATCH --job-name=run_HD
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=801-1000%10 # FIX
#SBATCH --nice=10000
#SBATCH --output=log/run_HD_801-1000-run_%a.log # FIX
#SBATCH --error=log/run_HD_801-1000-run_%a.log # FIX

# THIS IS A TEMPLATE THAT NEEDS TO BE EDITTED FOR RUNNING HD ON ALL SAMPLES IN A JOB ARRAY WITH JOB BATCH WITH N TASKS

# Define paths
BAM_PATH="/data/target_id/dna/results/dnaseq_results/BAM"
TARGETS_PATH="/data/target_id/dna/results/dnaseq_results/CNV_VCF/cnvkit_results/imp_file_for_cnvkit/wgs-bin-size-1000bp/wgs_targets.bed"
ANTITARGETS_PATH="/data/target_id/dna/results/dnaseq_results/CNV_VCF/cnvkit_results/imp_file_for_cnvkit/wgs-bin-size-1000bp/wgs_antitargets.bed"
REFERENCE_CNN="/data/target_id/dna/results/dnaseq_results/CNV_VCF/cnvkit_results/imp_file_for_cnvkit/wgs-bin-size-1000bp/wgs_cnv_reference_200C.cnn"
OUTPUT_PATH="/data/target_id/dna/results/dnaseq_results/CNV_VCF/cnvkit_results/wgs-bin-size-1000bp_results"
BAM_LIST="/data/target_id/dna/results/dnaseq_results/CNV_VCF/sample_lists/all_samples_randomised.txt"
snp_vcf="/data/target_id/dna/results/dnaseq_results/full_cohort_VCFs/04_annotated/relabelled_annotated/merged_relabelled_vcfs_30_03_23.vcf"

# Read BAM filename from the list based on the SLURM array job index
BAM_FILE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BAM_LIST)
BAM_FILE="$BAM_PATH/$BAM_FILE_NAME"
BASE_NAME=$(basename $BAM_FILE .bam)

# Extract sample ID from BAM file name
sample_id=$(echo $BAM_FILE_NAME | sed 's/\.post\.sortNMMD\.bam//')

# Ensure output directory exists
mkdir -p $OUTPUT_PATH/$BASE_NAME

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

# Step 1: Coverage calculation
echo "Starting coverage calculations for $BASE_NAME..."
run_step "$OUTPUT_PATH/$BASE_NAME/$BASE_NAME.targetcoverage.cnn" \
    "Calculating target coverage..." \
    "Target coverage calculation complete." \
    "cnvkit.py coverage $BAM_FILE $TARGETS_PATH -o $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.targetcoverage.cnn"

run_step "$OUTPUT_PATH/$BASE_NAME/$BASE_NAME.antitargetcoverage.cnn" \
    "Calculating antitarget coverage..." \
    "Antitarget coverage calculation complete." \
    "cnvkit.py coverage $BAM_FILE $ANTITARGETS_PATH -o $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.antitargetcoverage.cnn"

# Step 2: Fixing the coverage
run_step "$OUTPUT_PATH/$BASE_NAME/$BASE_NAME.cnr" \
    "Fixing coverage..." \
    "Coverage fixed." \
    "cnvkit.py fix $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.targetcoverage.cnn $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.antitargetcoverage.cnn $REFERENCE_CNN -o $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.cnr"

# Step 3: Segmentation
run_step "$OUTPUT_PATH/$BASE_NAME/$BASE_NAME.cns" \
    "Segmenting..." \
    "Segmentation complete." \
    "cnvkit.py segment $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.cnr -o $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.cns"

# Step 4: Calling copy number alterations with SNPs
run_step "$OUTPUT_PATH/$BASE_NAME/${BASE_NAME}.snp.cns" \
    "Calling CNVs with SNPs..." \
    "VCF generation complete." \
    "cnvkit.py call $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.cns -v $snp_vcf -i $sample_id -o $OUTPUT_PATH/$BASE_NAME/${BASE_NAME}.snp.cns"

# Step 5: Constructing VCF
run_step "$OUTPUT_PATH/$BASE_NAME/${BASE_NAME}.snp.vcf" \
    "Exporting CNVs to VCF format..." \
    "VCF export complete." \
    "cnvkit.py export vcf $OUTPUT_PATH/$BASE_NAME/$BASE_NAME.snp.cns -o $OUTPUT_PATH/$BASE_NAME/${BASE_NAME}.snp.vcf"

