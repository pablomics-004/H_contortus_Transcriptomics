#!/usr/bin/env bash

# ============================================================================ #
#  RNA-Seq Read Trimming Pipeline using Cutadapt
# ---------------------------------------------------------------------------- #
#  Description:
#    Trims RNA-Seq reads (paired-end or single-end) using Cutadapt.
#    Supports nested directory structure: SRA/SRS_ID/*.fastq.gz
#    Outputs all trimmed files into a single output directory.
#
#  Usage:
#    ./trimming.sh [DATADIR] [OUTDIR] [NCORES] [QUALITY] [MINLEN]
#
#  Default parameters:
#    DATADIR  = data
#    OUTDIR   = results/TRIMMED
#    NCORES   = 4
#    QUALITY  = 20
#    MINLEN   = 25
#
#  Author: Ashley Yael Montiel Vargas
#  Date: 2025-OCT
# ============================================================================ #

set -euo pipefail
shopt -s nullglob

# ============================== #
# Configuration / Parameters
# ============================== #
DATADIR="${1:-data}"
OUTDIR="${2:-results/TRIMMED}"
NCORES="${3:-4}"
QUALITY="${4:-20}"
MINLEN="${5:-25}"

mkdir -p "$OUTDIR"
LOG_FILE="$OUTDIR/cutadapt_trimming.log"

exec > >(tee -a "$LOG_FILE") 2>&1

# ============================== #
# Logging helpers
# ============================== #

# -------------------------------------------------------------
# log
# -------------------------------------------------------------
# Description:
#   Print timestamped messages to stdout and log file.
# Parameters:
#   $* - Message to print
# Output:
#   Prints message with timestamp
log() {
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# -------------------------------------------------------------
# error_exit
# -------------------------------------------------------------
# Description:
#   Print an error message to stderr and exit the script with code 1.
# Parameters:
#   $* - Error message
# Output:
#   Script exits immediately
error_exit() {
    echo -e "[ERROR] $*" >&2
    exit 1
}

# ============================== #
# Input validation
# ============================== #

# -------------------------------------------------------------
# validate_inputs
# -------------------------------------------------------------
# Description:
#   Validate input and output directories exist or can be created.
#   Validate that at least one FASTQ file exists in the data directory.
# Parameters:
#   DATADIR, OUTDIR
# Output:
#   Logs status; exits if validation fails
validate_inputs() {
    [[ -d "$DATADIR" ]] || error_exit "Data directory not found: $DATADIR"
    mkdir -p "$OUTDIR" || error_exit "Cannot create output directory: $OUTDIR"

    log "Input directory: $DATADIR"
    log "Output directory: $OUTDIR"
    log "Cores: $NCORES | Quality cutoff: $QUALITY | Minimum length: $MINLEN"

    FILES=( "$DATADIR"/*/*.fastq.gz )
    [[ ${#FILES[@]} -gt 0 ]] || error_exit "No FASTQ files found in $DATADIR"
}

# ============================== #
# Collect R1 FASTQ files for trimming
# ============================== #

# -------------------------------------------------------------
# collect_r1_files
# -------------------------------------------------------------
# Description:
#   Collect FASTQ files representing R1 reads for paired-end samples.
#   Single-end reads will also be included.
# Parameters:
#   FILES array (global)
# Output:
#   ITERATOR array with R1 files to process
collect_r1_files() {
    ITERATOR=()
    for f in "${FILES[@]}"; do
        [[ "$f" == *_2.fastq.gz ]] && continue
        ITERATOR+=( "$f" )
    done

    [[ ${#ITERATOR[@]} -gt 0 ]] || error_exit "No R1 FASTQ files found for processing"
    log "Found ${#ITERATOR[@]} R1 files for trimming"
}

# ============================== #
# Trimming function
# ============================== #

# -------------------------------------------------------------
# trim_reads
# -------------------------------------------------------------
# Description:
#   Runs Cutadapt on each R1 file, detecting paired-end or single-end automatically.
#   Writes all trimmed files to the output directory.
# Parameters:
#   ITERATOR array (R1 files), NCORES, QUALITY, MINLEN, OUTDIR
# Output:
#   Trimmed FASTQ files in OUTDIR; logs progress and errors
trim_reads() {
    log "=== Starting Cutadapt trimming ==="
    for r1 in "${ITERATOR[@]}"; do
        base=$(basename "$r1")
        r2="${r1/_1.fastq.gz/_2.fastq.gz}"

        clean_r1="$OUTDIR/${base/_1.fastq.gz/_clean_1.fastq.gz}"
        clean_r2="$OUTDIR/${base/_1.fastq.gz/_clean_2.fastq.gz}"
        clean_se="$OUTDIR/${base/.fastq.gz/_clean.fastq.gz}"

        if [[ -s "$r2" ]]; then
            log "[PE] Trimming $r1 & $r2"
            cutadapt -q "$QUALITY" -m "$MINLEN" -j "$NCORES" -o "$clean_r1" -p "$clean_r2" "$r1" "$r2" \
                || error_exit "Cutadapt failed on paired-end files: $r1, $r2"
        else
            log "[SE] Trimming $r1"
            cutadapt -q "$QUALITY" -m "$MINLEN" -j "$NCORES" -o "$clean_se" "$r1" \
                || error_exit "Cutadapt failed on single-end file: $r1"
        fi
    done
    log "[OK] Cutadapt trimming finished successfully. Trimmed files are in: $OUTDIR"
}

# ============================== #
# Main workflow
# ============================== #
main() {
    validate_inputs
    collect_r1_files
    trim_reads
}

main "$@"

