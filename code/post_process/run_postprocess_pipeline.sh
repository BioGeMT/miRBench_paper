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
    cat <<USAGE
Usage:
  Cohort mode:
    $0 --mode cohort -i input_dir -o output_dir [-p phyloP_file -c phastCons_file] [-r conservation_dir]

  Single mode:
    $0 --mode single -f input_file -o output_dir [-p phyloP_file -c phastCons_file] [-r conservation_dir]

Notes:
  - Provide both -p and -c, or omit both to use cached/downloaded conservation files.
  - Cohort mode expects the three canonical input datasets in input_dir.
  - Single mode processes exactly one dataset file.
USAGE
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

resolve_abs_path() {
    local path="$1"
    echo "$(cd "$(dirname "$path")" && pwd)/$(basename "$path")"
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

prepare_conservation_files() {
    conservation_dir=${conservation_dir:-$output_dir/reference_data/conservation}
    mkdir -p "$conservation_dir"
    conservation_dir=$(cd "$conservation_dir" && pwd)

    if [ -n "${phyloP_file:-}" ] && [ -n "${phastCons_file:-}" ]; then
        require_file "$phyloP_file"
        require_file "$phastCons_file"
        phyloP_file=$(resolve_abs_path "$phyloP_file")
        phastCons_file=$(resolve_abs_path "$phastCons_file")
        return
    fi

    phyloP_file="$conservation_dir/$PHYLOP_FILENAME"
    phastCons_file="$conservation_dir/$PHASTCONS_FILENAME"

    echo "Conservation files not provided explicitly. Using cache/download directory: $conservation_dir"
    ensure_conservation_file "$phyloP_file" "$PHYLOP_URL"
    ensure_conservation_file "$phastCons_file" "$PHASTCONS_URL"
}

run_postprocess0() {
    local input_dir="$1"
    local output_dir="$2"
    local intermediate_dir="$3"

    run_step "Step 0: filter and deduplicate" \
        bash "$SCRIPT_DIR/postprocess_0_filter_and_deduplicate.sh" \
        -i "$input_dir" \
        -o "$output_dir" \
        -n "$intermediate_dir"
}

run_postprocess1() {
    local input_dir="$1"
    local output_dir="$2"

    run_step "Step 1: genomic region annotation" \
        bash "$SCRIPT_DIR/postprocess_1_genomic_region_annotation.sh" \
        -i "$input_dir" \
        -o "$output_dir"
}

run_postprocess2() {
    local input_file="$1"
    local input_file_2="$2"
    local input_file_3="$3"
    local output_dir="$4"
    local intermediate_dir="$5"

    run_step "Step 2: exclude miRNA families" \
        bash "$SCRIPT_DIR/postprocess_2_exclude_mirna_families.sh" \
        -i "$input_file" \
        -j "$input_file_2" \
        -k "$input_file_3" \
        -o "$output_dir" \
        -n "$intermediate_dir"
}

run_postprocess3() {
    local input_dir="$1"
    local output_dir="$2"
    local intermediate_dir="$3"

    run_step "Step 3: make negatives" \
        bash "$SCRIPT_DIR/postprocess_3_make_negatives.sh" \
        -i "$input_dir" \
        -o "$output_dir" \
        -n "$intermediate_dir"
}

run_postprocess4() {
    local input_dir="$1"
    local output_dir="$2"

    run_step "Step 4: train/test splits" \
        bash "$SCRIPT_DIR/postprocess_4_train_test_splits.sh" \
        -i "$input_dir" \
        -o "$output_dir"
}

run_postprocess5() {
    local input_dir="$1"
    local output_dir="$2"

    run_step "Step 5: drop test column" \
        bash "$SCRIPT_DIR/postprocess_5_drop_test_col.sh" \
        -i "$input_dir" \
        -o "$output_dir"
}

run_postprocess6() {
    local input_dir="$1"
    local output_dir="$2"

    run_step "Step 6: add conservation" \
        bash "$SCRIPT_DIR/postprocess_6_add_conservation.sh" \
        -i "$input_dir" \
        -o "$output_dir" \
        -p "$phyloP_file" \
        -c "$phastCons_file"
}

publish_gzipped_output() {
    local source_file="$1"
    local target_file="$2"

    require_file "$source_file"
    mkdir -p "$(dirname "$target_file")"
    echo "Publishing $(basename "$source_file") -> $(basename "$target_file")"
    gzip -c "$source_file" > "$target_file"
    if [ ! -f "$target_file" ]; then
        echo "Error: Failed to create $target_file"
        exit 1
    fi
}

publish_cohort_outputs() {
    local step6_dir="$1"
    local release_dir="$output_dir/zenodo_release"

    local manakov_leftout_src="$step6_dir/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.excluded.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.conservation.tsv"
    local manakov_train_src="$step6_dir/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.remaining.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.conservation.tsv"
    local manakov_test_src="$step6_dir/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.remaining.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.conservation.tsv"
    local klimentova_test_src="$step6_dir/${KLIMENTOVA_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.conservation.tsv"
    local hejret_train_src="$step6_dir/${HEJRET_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.conservation.tsv"
    local hejret_test_src="$step6_dir/${HEJRET_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.conservation.tsv"

    mkdir -p "$release_dir"

    publish_gzipped_output "$manakov_leftout_src" "$release_dir/AGO2_eCLIP_Manakov2022_leftout.tsv.gz"
    publish_gzipped_output "$manakov_test_src" "$release_dir/AGO2_eCLIP_Manakov2022_test.tsv.gz"
    publish_gzipped_output "$manakov_train_src" "$release_dir/AGO2_eCLIP_Manakov2022_train.tsv.gz"
    publish_gzipped_output "$klimentova_test_src" "$release_dir/AGO2_eCLIP_Klimentova2022_test.tsv.gz"
    publish_gzipped_output "$hejret_test_src" "$release_dir/AGO2_CLASH_Hejret2023_test.tsv.gz"
    publish_gzipped_output "$hejret_train_src" "$release_dir/AGO2_CLASH_Hejret2023_train.tsv.gz"

    echo
    echo "Published Zenodo-ready files to: $release_dir"
}

run_cohort_mode() {
    local step0_dir="$output_dir/step0_filter_and_deduplicate"
    local step0_intermediate_dir="$step0_dir/intermediate"
    local step1_dir="$output_dir/step1_genomic_region_annotation"
    local step2_dir="$output_dir/step2_exclude_mirna_families"
    local step2_intermediate_dir="$step2_dir/intermediate"
    local step3_input_dir="$output_dir/step3_make_negatives/input"
    local step3_dir="$output_dir/step3_make_negatives/output"
    local step3_intermediate_dir="$step3_dir/intermediate"
    local step4_input_dir="$output_dir/step4_train_test_splits/input"
    local step4_dir="$output_dir/step4_train_test_splits/output"
    local step5_input_dir="$output_dir/step5_drop_test_col/input"
    local step5_dir="$output_dir/step5_drop_test_col/output"
    local step6_dir="$output_dir/step6_add_conservation"

    local manakov_input="$input_dir/${MANAKOV_BASE}.tsv"
    local hejret_input="$input_dir/${HEJRET_BASE}.tsv"
    local klimentova_input="$input_dir/${KLIMENTOVA_BASE}.tsv"

    require_file "$manakov_input"
    require_file "$hejret_input"
    require_file "$klimentova_input"

    mkdir -p \
        "$step0_dir" \
        "$step0_intermediate_dir" \
        "$step1_dir" \
        "$step2_dir" \
        "$step2_intermediate_dir" \
        "$step3_input_dir" \
        "$step3_dir" \
        "$step3_intermediate_dir" \
        "$step4_input_dir" \
        "$step4_dir" \
        "$step5_input_dir" \
        "$step5_dir" \
        "$step6_dir"

    run_postprocess0 "$input_dir" "$step0_dir" "$step0_intermediate_dir"

    local manakov_step0="$step0_dir/${MANAKOV_BASE}.filtered.deduplicated.tsv"
    local hejret_step0="$step0_dir/${HEJRET_BASE}.filtered.deduplicated.tsv"
    local klimentova_step0="$step0_dir/${KLIMENTOVA_BASE}.filtered.deduplicated.tsv"

    require_file "$manakov_step0"
    require_file "$hejret_step0"
    require_file "$klimentova_step0"

    run_postprocess1 "$step0_dir" "$step1_dir"

    local manakov_step1="$step1_dir/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.tsv"
    local hejret_step1="$step1_dir/${HEJRET_BASE}.filtered.deduplicated.annotated.filtered.tsv"
    local klimentova_step1="$step1_dir/${KLIMENTOVA_BASE}.filtered.deduplicated.annotated.filtered.tsv"

    require_file "$manakov_step1"
    require_file "$hejret_step1"
    require_file "$klimentova_step1"

    run_postprocess2 "$manakov_step1" "$hejret_step1" "$klimentova_step1" "$step2_dir" "$step2_intermediate_dir"

    local manakov_excluded="$step2_dir/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.excluded.tsv"
    local manakov_remaining="$step2_dir/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.remaining.tsv"

    require_file "$manakov_excluded"
    require_file "$manakov_remaining"

    link_into_dir "$manakov_remaining" "$step3_input_dir"
    link_into_dir "$manakov_excluded" "$step3_input_dir"
    link_into_dir "$hejret_step1" "$step3_input_dir"
    link_into_dir "$klimentova_step1" "$step3_input_dir"

    run_postprocess3 "$step3_input_dir" "$step3_dir" "$step3_intermediate_dir"

    local manakov_neg="$step3_dir/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.remaining.gene_clusters_added.mirfam_sorted.negatives.tsv"
    local hejret_neg="$step3_dir/${HEJRET_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.tsv"
    local klimentova_neg="$step3_dir/${KLIMENTOVA_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.tsv"
    local manakov_excluded_neg="$step3_dir/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.excluded.gene_clusters_added.mirfam_sorted.negatives.tsv"

    require_file "$manakov_neg"
    require_file "$hejret_neg"
    require_file "$klimentova_neg"
    require_file "$manakov_excluded_neg"

    link_into_dir "$manakov_neg" "$step4_input_dir"
    link_into_dir "$hejret_neg" "$step4_input_dir"

    run_postprocess4 "$step4_input_dir" "$step4_dir"

    local manakov_train="$step4_dir/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.remaining.gene_clusters_added.mirfam_sorted.negatives.train.tsv"
    local manakov_test="$step4_dir/${MANAKOV_BASE}.filtered.deduplicated.annotated.filtered.remaining.gene_clusters_added.mirfam_sorted.negatives.test.tsv"
    local hejret_train="$step4_dir/${HEJRET_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.train.tsv"
    local hejret_test="$step4_dir/${HEJRET_BASE}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.test.tsv"

    require_file "$manakov_train"
    require_file "$manakov_test"
    require_file "$hejret_train"
    require_file "$hejret_test"

    link_into_dir "$manakov_train" "$step5_input_dir"
    link_into_dir "$manakov_test" "$step5_input_dir"
    link_into_dir "$hejret_train" "$step5_input_dir"
    link_into_dir "$hejret_test" "$step5_input_dir"
    link_into_dir "$klimentova_neg" "$step5_input_dir"
    link_into_dir "$manakov_excluded_neg" "$step5_input_dir"

    run_postprocess5 "$step5_input_dir" "$step5_dir"
    run_postprocess6 "$step5_dir" "$step6_dir"
    publish_cohort_outputs "$step6_dir"

    echo
    echo "Cohort mode completed successfully."
    echo "Final outputs are in: $step6_dir"
    echo "Zenodo-ready files are in: $output_dir/zenodo_release"
}

run_single_mode() {
    local step0_dir="$output_dir/step0_filter_and_deduplicate"
    local step0_input_dir="$output_dir/step0_filter_and_deduplicate/input"
    local step0_intermediate_dir="$step0_dir/intermediate"
    local step1_dir="$output_dir/step1_genomic_region_annotation"
    local step3_dir="$output_dir/step3_make_negatives"
    local step3_intermediate_dir="$step3_dir/intermediate"
    local step4_dir="$output_dir/step4_train_test_splits"
    local step5_dir="$output_dir/step5_drop_test_col"
    local step6_dir="$output_dir/step6_add_conservation"

    require_file "$input_file"
    local input_file_abs
    input_file_abs=$(resolve_abs_path "$input_file")
    local base_name
    base_name=$(basename "$input_file_abs" .tsv)

    mkdir -p \
        "$step0_input_dir" \
        "$step0_intermediate_dir" \
        "$step1_dir" \
        "$step3_dir" \
        "$step3_intermediate_dir" \
        "$step4_dir" \
        "$step5_dir" \
        "$step6_dir"

    link_into_dir "$input_file_abs" "$step0_input_dir"

    run_postprocess0 "$step0_input_dir" "$step0_dir" "$step0_intermediate_dir"

    local step0_file="$step0_dir/${base_name}.filtered.deduplicated.tsv"
    require_file "$step0_file"

    run_postprocess1 "$step0_dir" "$step1_dir"

    local step1_file="$step1_dir/${base_name}.filtered.deduplicated.annotated.filtered.tsv"
    require_file "$step1_file"

    run_postprocess3 "$step1_dir" "$step3_dir" "$step3_intermediate_dir"

    local negatives_file="$step3_dir/${base_name}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.tsv"
    require_file "$negatives_file"

    run_postprocess4 "$step3_dir" "$step4_dir"

    local train_file="$step4_dir/${base_name}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.train.tsv"
    local test_file="$step4_dir/${base_name}.filtered.deduplicated.annotated.filtered.gene_clusters_added.mirfam_sorted.negatives.test.tsv"
    require_file "$train_file"
    require_file "$test_file"

    run_postprocess5 "$step4_dir" "$step5_dir"
    run_postprocess6 "$step5_dir" "$step6_dir"

    echo
    echo "Single mode completed successfully."
    echo "Final outputs are in: $step6_dir"
}

mode=""
input_dir=""
input_file=""
output_dir=""
phyloP_file=""
phastCons_file=""
conservation_dir=""

while [ $# -gt 0 ]; do
    case "$1" in
        --mode|-m)
            mode="$2"
            shift 2
            ;;
        --input-dir|-i)
            input_dir="$2"
            shift 2
            ;;
        --input-file|-f)
            input_file="$2"
            shift 2
            ;;
        --output-dir|-o)
            output_dir="$2"
            shift 2
            ;;
        --phyloP|-p)
            phyloP_file="$2"
            shift 2
            ;;
        --phastCons|-c)
            phastCons_file="$2"
            shift 2
            ;;
        --conservation-dir|-r)
            conservation_dir="$2"
            shift 2
            ;;
        --help|-h)
            usage
            ;;
        *)
            echo "Error: Unknown argument: $1"
            usage
            ;;
    esac
done

if [ -z "$mode" ] || [ -z "$output_dir" ]; then
    usage
fi

if [ "$mode" != "cohort" ] && [ "$mode" != "single" ]; then
    echo "Error: --mode must be 'cohort' or 'single'."
    usage
fi

if { [ -n "$phyloP_file" ] && [ -z "$phastCons_file" ]; } || { [ -z "$phyloP_file" ] && [ -n "$phastCons_file" ]; }; then
    echo "Error: Provide both -p and -c, or omit both."
    usage
fi

mkdir -p "$output_dir"
output_dir=$(cd "$output_dir" && pwd)

log_file="$output_dir/run_postprocess_pipeline.log"
exec > >(tee -a "$log_file") 2>&1

prepare_conservation_files

case "$mode" in
    cohort)
        if [ -z "$input_dir" ] || [ -n "$input_file" ]; then
            echo "Error: Cohort mode requires -i input_dir and does not use -f input_file."
            usage
        fi
        require_dir "$input_dir"
        input_dir=$(cd "$input_dir" && pwd)
        run_cohort_mode
        ;;
    single)
        if [ -z "$input_file" ] || [ -n "$input_dir" ]; then
            echo "Error: Single mode requires -f input_file and does not use -i input_dir."
            usage
        fi
        run_single_mode
        ;;
esac
