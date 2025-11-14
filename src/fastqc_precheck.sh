#!/usr/bin/env bash

# ============================================================================ #
#  RNA-Seq Preliminary Quality Check using FastQC
# ---------------------------------------------------------------------------- #
#  Description:
#    Runs FastQC on RNA-Seq FASTQ files to assess quality before trimming.
#    Supports nested directory structure: DATADIR/SRS_ID/*.fastq.gz
#    All FastQC outputs (HTML + ZIP) are written into a single output directory.
#
#  Usage:
#    ./fastqc_precheck.sh [DATADIR] [OUTDIR] [NCORES]
#
#  Parameters:
#    DATADIR  - Directory with raw FASTQ files (default: data)
#    OUTDIR   - Directory to store FastQC reports (default: results/FASTQC)
#    NCORES   - Number of CPU cores for FastQC (default: 4)
#
#  Author: Ashley Yael Montiel Vargas
#  Date: 2025-NOV
# ============================================================================ #


set -euo pipefail
shopt -s nullglob

# ============================== #
# Configuration / Parameters
# ============================== #
DATADIR="${1:-data}"
OUTDIR="${2:-results/FASTQC}"
NCORES="${3:-4}"

mkdir -p "$OUTDIR"
LOG_FILE="$OUTDIR/fastqc_precheck.log"

exec > >(tee -a "$LOG_FILE") 2>&1

# ============================== #
# Logging helpers
# ============================== #
log() {
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

error_exit() {
    echo -e "[ERROR] $*" >&2
    exit 1
}

# ============================== #
# Validate inputs
# ============================== #
validate_inputs() {
    [[ -d "$DATADIR" ]] || error_exit "Data directory not found: $DATADIR"
    mkdir -p "$OUTDIR" || error_exit "Cannot create output directory: $OUTDIR"
    log "Input directory: $DATADIR"
    log "Output directory: $OUTDIR"
    log "Cores: $NCORES"

    FILES=( "$DATADIR"/*/*.fastq.gz )
    [[ ${#FILES[@]} -gt 0 ]] || error_exit "No FASTQ files found in $DATADIR"
    log "Found ${#FILES[@]} FASTQ files for FastQC"
}

# ============================== #
# Collect files
# ============================== #
collect_fastq_files() {
    ITERATOR=()
    for f in "${FILES[@]}"; do
        ITERATOR+=( "$f" )
    done
    [[ ${#ITERATOR[@]} -gt 0 ]] || error_exit "No FASTQ files to process"
    log "Prepared ${#ITERATOR[@]} files for FastQC"
}

# ============================== #
# Run FastQC
# ============================== #
run_fastqc() {
    log "=== Starting FastQC analysis ==="
    for fq in "${ITERATOR[@]}"; do
        sample=$(basename "$fq")
        log "Processing $sample"
        fastqc -o "$OUTDIR" -t "$NCORES" "$fq" \
            || error_exit "FastQC failed for file: $fq"
    done
    log "[OK] FastQC analysis completed. Reports are in: $OUTDIR"
}

# ============================== #
# Main workflow
# ============================== #
main() {
    validate_inputs
    collect_fastq_files
    run_fastqc
}

main "$@"
