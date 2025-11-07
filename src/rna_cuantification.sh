#!/usr/bin/env bash

# ============================================================================ #
#  RNA-Seq Quantification Pipeline: Salmon or Kallisto
# ---------------------------------------------------------------------------- #
#  Description:
#    Modular script for transcript quantification using Salmon or Kallisto.
#    It builds an index (if missing), quantifies samples (paired/single-end),
#    and produces a combined count matrix.
#
#  Usage:
#    ./rna_quantification.sh [TOOL] [TRIM_DIR] [TRANSCRIPTOME_FASTA] [OUTDIR] [NCORES]
#
#  Example:
#    ./rna_quantification.sh salmon  results/TRIMMED data/Trinity.fasta results/SALMON   8
#    ./rna_quantification.sh kallisto results/TRIMMED data/Trinity.fasta results/KALLISTO 8
#
#  Environment paths:
#    Modify these if tools are installed in a specific environment.
#    Example (Bioconda):
#      SALMON_EXE="/export/apps/bioconda/envs/salmon/bin/salmon"
#      KALLISTO_EXE="/export/apps/bioconda/envs/kallisto/bin/kallisto"
#
#  Author: Ashley Yael Montiel Vargas
#  Date: 2025-OCT
# ============================================================================ #

set -euo pipefail
shopt -s nullglob

# ----------------------------- #
#  Configuration: Tool paths
# ----------------------------- #
SALMON_EXE=${SALMON_EXE:-/export/apps/bioconda/envs/salmon/bin/salmon}
KALLISTO_EXE=${KALLISTO_EXE:-/export/apps/bioconda/envs/kallisto/bin/kallisto}

# ----------------------------- #
#  Default arguments
# ----------------------------- #
TOOL="${1:-salmon}"                              # salmon | kallisto
TRIM_DIR="${2:-results/TRIMMED}"
TRANSCRIPTOME_FASTA="${3:-data/Trinity.fasta}"
OUTDIR="${4:-results/${TOOL^^}}"
NCORES="${5:-8}"

mkdir -p "$OUTDIR"
LOG_FILE="$OUTDIR/rna_cuantification.log"

exec > >(tee -a "$LOG_FILE") 2>&1

# ----------------------------- #
#  Logging helpers
# ----------------------------- #
log() {
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

error_exit() {
    echo -e "[ERROR] $*" >&2
    exit 1
}

# ----------------------------- #
#  Validation functions
# ----------------------------- #

check_inputs() {
    log "Validating inputs..."
    [[ -d "$TRIM_DIR" ]] || error_exit "Trimmed reads directory not found: $TRIM_DIR"
    [[ -f "$TRANSCRIPTOME_FASTA" ]] || error_exit "Transcriptome FASTA not found: $TRANSCRIPTOME_FASTA"
    [[ "$NCORES" =~ ^[0-9]+$ ]] || error_exit "Invalid number of cores: $NCORES"
}

check_executables() {
    log "Checking required executables..."
    if [[ "$TOOL" == "salmon" ]]; then
        [[ -x "$SALMON_EXE" ]] || error_exit "Salmon executable not found: $SALMON_EXE"
    elif [[ "$TOOL" == "kallisto" ]]; then
        [[ -x "$KALLISTO_EXE" ]] || error_exit "Kallisto executable not found: $KALLISTO_EXE"
    else
        error_exit "Invalid tool specified: $TOOL (must be 'salmon' or 'kallisto')"
    fi
}

check_fastq_files() {
    log "Checking FASTQ files..."
    files=("$TRIM_DIR"/*.fastq "$TRIM_DIR"/*.fastq.gz)
    if [[ ${#files[@]} -eq 0 ]]; then
        error_exit "No FASTQ files found in $TRIM_DIR"
    fi
    log "Found ${#files[@]} FASTQ files."
}

# ----------------------------- #
#  Detect paired and single-end
# ----------------------------- #
detect_pairs() {
    log "Detecting paired- and single-end samples..."
    declare -gA PE_R1=()
    declare -gA PE_R2=()
    declare -gA SE=()

    for f in "$TRIM_DIR"/*.fastq "$TRIM_DIR"/*.fastq.gz; do
        [[ -f "$f" ]] || continue
        base=$(basename "$f")

        # Common patterns for read 1 and read 2
        if [[ "$base" =~ (_R?1[^0-9]*\.f(ast)?q(\.gz)?)$ ]]; then
            sample="${base%%_R*}"
            r2_candidate="${f/_R1/_R2}"
            r2_candidate="${r2_candidate/_1/_2}"
            if [[ -f "$r2_candidate" ]]; then
                PE_R1["$sample"]="$f"
                PE_R2["$sample"]="$r2_candidate"
            else
                SE["$sample"]="$f"
            fi
        elif [[ "$base" =~ (_R?2[^0-9]*\.f(ast)?q(\.gz)?)$ ]]; then
            continue  # already paired
        else
            # No _1 or _2 pattern — treat as single-end
            sample="${base%%.*}"
            SE["$sample"]="$f"
        fi
    done

    log "Detected ${#PE_R1[@]} paired-end and ${#SE[@]} single-end samples."
}

# ----------------------------- #
#  Index creation
# ----------------------------- #
build_index_salmon() {
    SALMON_INDEX="$OUTDIR/salmon_index"
    if [[ ! -d "$SALMON_INDEX" ]]; then
        log "Building Salmon index..."
        "$SALMON_EXE" index -t "$TRANSCRIPTOME_FASTA" -i "$SALMON_INDEX" -k 31
    else
        log "Salmon index already exists: $SALMON_INDEX"
    fi
}

build_index_kallisto() {
    KALLISTO_INDEX="$OUTDIR/transcriptome.idx"
    if [[ ! -f "$KALLISTO_INDEX" ]]; then
        log "Building Kallisto index..."
        "$KALLISTO_EXE" index -i "$KALLISTO_INDEX" "$TRANSCRIPTOME_FASTA"
    else
        log "Kallisto index already exists: $KALLISTO_INDEX"
    fi
}

# ----------------------------- #
#  Quantification
# ----------------------------- #
quantify_salmon() {
    SALMON_INDEX="$OUTDIR/salmon_index"
    log "Running Salmon quantification..."

    for s in "${!PE_R1[@]}"; do
        out="$OUTDIR/$s"; mkdir -p "$out"
        log "[PE] Quantifying $s"
        "$SALMON_EXE" quant -i "$SALMON_INDEX" -l A \
            -1 "${PE_R1[$s]}" -2 "${PE_R2[$s]}" \
            -p "$NCORES" -o "$out"
    done

    for s in "${!SE[@]}"; do
        out="$OUTDIR/$s"; mkdir -p "$out"
        log "[SE] Quantifying $s"
        "$SALMON_EXE" quant -i "$SALMON_INDEX" -l A \
            -r "${SE[$s]}" \
            -p "$NCORES" -o "$out"
    done
}

quantify_kallisto() {
    KALLISTO_INDEX="$OUTDIR/transcriptome.idx"
    log "Running Kallisto quantification..."

    for s in "${!PE_R1[@]}"; do
        out="$OUTDIR/$s"; mkdir -p "$out"
        log "[PE] Quantifying $s"
        "$KALLISTO_EXE" quant -i "$KALLISTO_INDEX" -o "$out" -t "$NCORES" \
            "${PE_R1[$s]}" "${PE_R2[$s]}"
    done

    for s in "${!SE[@]}"; do
        out="$OUTDIR/$s"; mkdir -p "$out"
        log "[SE] Quantifying $s"
        "$KALLISTO_EXE" quant -i "$KALLISTO_INDEX" -o "$out" -t "$NCORES" \
            --single -l 200 -s 50 "${SE[$s]}"
    done
}

generate_matrix() {
    log "Generating count matrix..."
    COUNT_MATRIX="$OUTDIR/counts_matrix_${TOOL}.tsv"

    tmpdir=$(mktemp -d)
    samples=()

    # Select correct file pattern
    pattern="quant.sf"
    [[ "$TOOL" == "kallisto" ]] && pattern="abundance.tsv"

    # Extract sample names and headers
    for f in "$OUTDIR"/*/"$pattern"; do
        [[ -f "$f" ]] || continue
        sample=$(basename "$(dirname "$f")")
        samples+=("$sample")
        # Column with counts: 5th in Salmon, 4th in Kallisto
        col=$([[ "$TOOL" == "salmon" ]] && echo 5 || echo 4)
        awk -v c="$col" 'NR>1 {print $1"\t"$c}' "$f" > "$tmpdir/$sample.tsv"
    done

    [[ ${#samples[@]} -gt 0 ]] || error_exit "No quantification files found."

    # Merge all TSVs by transcript_id
    cp "$tmpdir/${samples[0]}.tsv" "$tmpdir/merged.tsv"
    for i in $(seq 1 $((${#samples[@]} - 1))); do
        join -t $'\t' "$tmpdir/merged.tsv" "$tmpdir/${samples[i]}.tsv" > "$tmpdir/tmp.tsv"
        mv "$tmpdir/tmp.tsv" "$tmpdir/merged.tsv"
    done

    # Write header and final matrix
    echo -e "transcript_id\t${samples[*]}" | tr ' ' '\t' > "$COUNT_MATRIX"
    cat "$tmpdir/merged.tsv" >> "$COUNT_MATRIX"
    rm -r "$tmpdir"

    log "Count matrix generated at: $COUNT_MATRIX"
}

# ----------------------------- #
#  Main workflow
# ----------------------------- #
main() {
    log "=== Starting RNA-Seq Quantification Pipeline ($TOOL) ==="
    check_inputs
    check_executables
    check_fastq_files
    detect_pairs

    if [[ "$TOOL" == "salmon" ]]; then
        build_index_salmon
        quantify_salmon
    else
        build_index_kallisto
        quantify_kallisto
    fi

    generate_matrix
    log "=== Pipeline completed successfully! Results in: $OUTDIR ==="
}

main "$@"
