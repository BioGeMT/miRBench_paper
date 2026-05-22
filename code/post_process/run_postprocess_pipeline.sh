#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=run_postprocess_pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

MANAKOV_BASE="AGO2_eCLIP_Manakov22_full_dataset"
HEJRET_BASE="AGO2_CLASH_Hejret2023_full_dataset"
KLIMENTOVA_BASE="AGO2_eCLIP_Klimentova22_full_dataset"
PHYLOP_FILENAME="hg38.phyloP100way.bw"
PHASTCONS_FILENAME="hg38.phastCons100way.bw"
PHYLOP_URL="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/${PHYLOP_FILENAME}"
PHASTCONS_URL="https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/${PHASTCONS_FILENAME}"

usage() {
    echo "Usage: $0 -i input_dir -o output_dir [-p phyloP_file -c phastCons_file] [-r conservation_dir]"
    echo "Provide both -p and -c, or omit both to use cached/downloaded files."
    exit 1
}

require_file() {
    local path="$1"
    if [ ! -f "$path" ]; then
        echo "Error: Required file not found: $path"
        exit 1
    fi
}

require_dir() {
    local path="$1"
    if [ ! -d "$path" ]; then
        echo "Error: Required directory not found: $path"
        exit 1
    fi
}

link_into_dir() {
    local source_path="$1"
    local target_dir="$2"
    require_file "$source_path"
    mkdir -p "$target_dir"
    ln -sfn "$source_path" "$target_dir/$(basename "$source_path")"
}

run_step() {
    local description="$1"
    shift
    echo
    echo "=== $description ==="
    "$@"
}

resolve_abs_path() {
    local path="$1"
    echo "$(cd "$(dirname "$path")" && pwd)/$(basename "$path")"
}

ensure_conservation_file() {
    local target_file="$1"
    local download_url="$2"

    if [ -f "$target_file" ]; then
        echo "Using cached conservation file: $target_file"
        return
    fi

    if ! command -v wget >/dev/null 2>&1; then
        echo "Error: wget is required to download conservation files automatically."
        exit 1
    fi

    mkdir -p "$(dirname "$target_file")"
    echo "Downloading $(basename "$target_file") ..."
    wget -O "$target_file" "$download_url"

    if [ ! -f "$target_file" ]; then
        echo "Error: Download failed for $target_file"
        exit 1
    fi
}

while getopts i:o:p:c:r: flag; do
    case "${flag}" in
        i) input_dir=${OPTARG} ;;
        o) output_dir=${OPTARG} ;;
        p) phyloP_file=${OPTARG} ;;
        c) phastCons_file=${OPTARG} ;;
        r) conservation_dir=${OPTARG} ;;
        *) usage ;;
    esac
done

if [ -z "${input_dir:-}" ] || [ -z "${output_dir:-}" ]; then
    usage
fi

if { [ -n "${phyloP_file:-}" ] && [ -z "${phastCons_file:-}" ]; } || { [ -z "${phyloP_file:-}" ] && [ -n "${phastCons_file:-}" ]; }; then
    echo "Error: Provide both -p and -c, or omit both."
    usage
fi

require_dir "$input_dir"

mkdir -p "$output_dir"

input_dir=$(cd "$input_dir" && pwd)
output_dir=$(cd "$output_dir" && pwd)
conservation_dir=${conservation_dir:-$output_dir/reference_data/conservation}
mkdir -p "$conservation_dir"
conservation_dir=$(cd "$conservation_dir" && pwd)

if [ -n "${phyloP_file:-}" ] && [ -n "${phastCons_file:-}" ]; then
    require_file "$phyloP_file"
    require_file "$phastCons_file"
    phyloP_file=$(resolve_abs_path "$phyloP_file")
    phastCons_file=$(resolve_abs_path "$phastCons_file")
else
    phyloP_file="$conservation_dir/$PHYLOP_FILENAME"
    phastCons_file="$conservation_dir/$PHASTCONS_FILENAME"
fi

log_file="$output_dir/run_postprocess_pipeline.log"
exec > >(tee -a "$log_file") 2>&1

if [ ! -f "$phyloP_file" ] || [ ! -f "$phastCons_file" ]; then
    echo "Conservation files not provided explicitly. Using cache/download directory: $conservation_dir"
    ensure_conservation_file "$phyloP_file" "$PHYLOP_URL"
    ensure_conservation_file "$phastCons_file" "$PHASTCONS_URL"
fi

MANAKOV_INPUT="$input_dir/${MANAKOV_BASE}.tsv"
HEJRET_INPUT="$input_dir/${HEJRET_BASE}.tsv"
KLIMENTOVA_INPUT="$input_dir/${KLIMENTOVA_BASE}.tsv"

require_file "$MANAKOV_INPUT"
require_file "$HEJRET_INPUT"
require_file "$KLIMENTOVA_INPUT"

STEP0_DIR="$output_dir/step0_filter_and_deduplicate"
STEP0_INTERMEDIATE_DIR="$STEP0_DIR/intermediate"
STEP1_DIR="$output_dir/step1_genomic_region_annotation"
STEP2_DIR="$output_dir/step2_exclude_mirna_families"
STEP2_INTERMEDIATE_DIR="$STEP2_DIR/intermediate"
STEP3_INPUT_DIR="$output_dir/step3_make_negatives/input"
STEP3_DIR="$output_dir/step3_make_negatives/output"
STEP3_INTERMEDIATE_DIR="$STEP3_DIR/intermediate"
STEP4_INPUT_DIR="$output_dir/step4_train_test_splits/input"
STEP4_DIR="$output_dir/step4_train_test_splits/output"
STEP5_INPUT_DIR="$output_dir/step5_drop_test_col/input"
STEP5_DIR="$output_dir/step5_drop_test_col/output"
STEP6_DIR="$output_dir/step6_add_conservation"

mkdir -p \
    "$STEP0_DIR" \
    "$STEP0_INTERMEDIATE_DIR" \
    "$STEP1_DIR" \
    "$STEP2_DIR" \
    "$STEP2_INTERMEDIATE_DIR" \
    "$STEP3_INPUT_DIR" \
    "$STEP3_DIR" \
    "$STEP3_INTERMEDIATE_DIR" \
    "$STEP4_INPUT_DIR" \
    "$STEP4_DIR" \
    "$STEP5_INPUT_DIR" \
    "$STEP5_DIR" \
    "$STEP6_DIR"

run_step "Step 0: filter and deduplicate" \
    bash "$SCRIPT_DIR/postprocess_0_filter_and_deduplicate.sh" \
    -i "$input_dir" \
    -o "$STEP0_DIR" \
    -n "$STEP0_INTERMEDIATE_DIR"

MANAKOV_STEP0="$STEP0_DIR/${MANAKOV_BASE}.filtered.deduplicated.tsv"
HEJRET_STEP0="$STEP0_DIR/${HEJRET_BASE}.filtered.deduplicated.tsv"
KLIMENTOVA_STEP0="$STEP0_DIR/${KLIMENTOVA_BASE}.filtered.deduplicated.tsv"

require_file "$MANAKOV_STEP0"
require_file "$HEJRET_STEP0"
require_file "$KLIMENTOVA_STEP0"

run_step "Step 1: genomic region annotation" \
    bash "$SCRIPT_DIR/postprocess_1_genomic_region_annotation.sh" \
    -i "$STEP0_DIR" \
    -o "$STEP1_DIR"

MANAKOV_STEP1="$STEP1_DIR/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.tsv"
HEJRET_STEP1="$STEP1_DIR/${HEJRET_BASE}.filtered.deduplicated.annotated.filtered.tsv"
KLIMENTOVA_STEP1="$STEP1_DIR/${KLIMENTOVA_BASE}.filtered.deduplicated.annotated.filtered.tsv"

require_file "$MANAKOV_STEP1"
require_file "$HEJRET_STEP1"
require_file "$KLIMENTOVA_STEP1"

run_step "Step 2: exclude miRNA families" \
    bash "$SCRIPT_DIR/postprocess_2_exclude_mirna_families.sh" \
    -i "$MANAKOV_STEP1" \
    -j "$HEJRET_STEP1" \
    -k "$KLIMENTOVA_STEP1" \
    -o "$STEP2_DIR" \
    -n "$STEP2_INTERMEDIATE_DIR"

MANAKOV_EXCLUDED="$STEP2_DIR/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.excluded.tsv"
MANAKOV_REMAINING="$STEP2_DIR/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.remaining.tsv"

require_file "$MANAKOV_EXCLUDED"
require_file "$MANAKOV_REMAINING"

link_into_dir "$MANAKOV_REMAINING" "$STEP3_INPUT_DIR"
link_into_dir "$MANAKOV_EXCLUDED" "$STEP3_INPUT_DIR"
link_into_dir "$HEJRET_STEP1" "$STEP3_INPUT_DIR"
link_into_dir "$KLIMENTOVA_STEP1" "$STEP3_INPUT_DIR"

run_step "Step 3: make negatives" \
    bash "$SCRIPT_DIR/postprocess_3_make_negatives.sh" \
    -i "$STEP3_INPUT_DIR" \
    -o "$STEP3_DIR" \
    -n "$STEP3_INTERMEDIATE_DIR"

MANAKOV_NEG="$STEP3_DIR/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.remaining.gene_clusters_added.mirfam_sorted.negatives.tsv"
HEJRET_NEG="$STEP3_DIR/${HEJRET_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.tsv"
KLIMENTOVA_NEG="$STEP3_DIR/${KLIMENTOVA_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.tsv"
MANAKOV_EXCLUDED_NEG="$STEP3_DIR/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.excluded.gene_clusters_added.mirfam_sorted.negatives.tsv"

require_file "$MANAKOV_NEG"
require_file "$HEJRET_NEG"
require_file "$KLIMENTOVA_NEG"
require_file "$MANAKOV_EXCLUDED_NEG"

link_into_dir "$MANAKOV_NEG" "$STEP4_INPUT_DIR"
link_into_dir "$HEJRET_NEG" "$STEP4_INPUT_DIR"

run_step "Step 4: train/test splits" \
    bash "$SCRIPT_DIR/postprocess_4_train_test_splits.sh" \
    -i "$STEP4_INPUT_DIR" \
    -o "$STEP4_DIR"

MANAKOV_TRAIN="$STEP4_DIR/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.remaining.gene_clusters_added.mirfam_sorted.negatives.train.tsv"
MANAKOV_TEST="$STEP4_DIR/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.remaining.gene_clusters_added.mirfam_sorted.negatives.test.tsv"
HEJRET_TRAIN="$STEP4_DIR/${HEJRET_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.train.tsv"
HEJRET_TEST="$STEP4_DIR/${HEJRET_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.test.tsv"

require_file "$MANAKOV_TRAIN"
require_file "$MANAKOV_TEST"
require_file "$HEJRET_TRAIN"
require_file "$HEJRET_TEST"

link_into_dir "$MANAKOV_TRAIN" "$STEP5_INPUT_DIR"
link_into_dir "$MANAKOV_TEST" "$STEP5_INPUT_DIR"
link_into_dir "$HEJRET_TRAIN" "$STEP5_INPUT_DIR"
link_into_dir "$HEJRET_TEST" "$STEP5_INPUT_DIR"
link_into_dir "$KLIMENTOVA_NEG" "$STEP5_INPUT_DIR"
link_into_dir "$MANAKOV_EXCLUDED_NEG" "$STEP5_INPUT_DIR"

run_step "Step 5: drop test column" \
    bash "$SCRIPT_DIR/postprocess_5_drop_test_col.sh" \
    -i "$STEP5_INPUT_DIR" \
    -o "$STEP5_DIR"

run_step "Step 6: add conservation" \
    bash "$SCRIPT_DIR/postprocess_6_add_conservation.sh" \
    -i "$STEP5_DIR" \
    -o "$STEP6_DIR" \
    -p "$phyloP_file" \
    -c "$phastCons_file"

echo
echo "run_postprocess_pipeline completed successfully."
echo "Final outputs are in: $STEP6_DIR"
