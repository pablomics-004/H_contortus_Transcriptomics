#!/usr/bin/env bash
log "=== [3/6] Alignment ==="

for c in "$TRIM_DIR"/*_clean*.fastq.gz; do
    base=$(basename "$c")
    sample="${base%%_*}"

    if [[ -f "$TRIM_DIR/${sample}_clean_2.fastq.gz" ]]; then
        bwa mem -t 8  "$BWA_DIR/BWA_INDEX/$BWA_PREFIX" \
            "$TRIM_DIR/${sample}_clean_1.fastq.gz" \
            "$TRIM_DIR/${sample}_clean_2.fastq.gz" \
            > "$ALIGN_DIR/${sample}.sam"
    else
        bwa mem -t 8 "$BWA_DIR/BWA_INDEX/$BWA_PREFIX" \
            "$TRIM_DIR/${sample}_clean.fastq.gz" \
            > "$ALIGN_DIR/${sample}.sam"
    fi

    wait_for_file "$ALIGN_DIR/${sample}.sam" "align $sample"
done
