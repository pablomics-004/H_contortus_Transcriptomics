#!/usr/bin/env bash
log "=== [1/6] Trimming ==="

declare -A samples
for f in "$DATADIR"/*.fastq.gz; do
    base=$(basename "$f")
    if [[ "$base" == *_1.fastq.gz ]]; then
        sample="${base%_1.fastq.gz}"
        samples[$sample]="PE"
    elif [[ "$base" != *_2.fastq.gz" ]]; then
        sample="${base%.fastq.gz}"
        samples[$sample]="SE"
    fi
done

for s in "${!samples[@]}"; do
    if [[ "${samples[$s]}" == "PE" ]]; then
        cutadapt -q 20 -m 25 \
            -o "$TRIM_DIR/${s}_clean_1.fastq.gz" \
            -p "$TRIM_DIR/${s}_clean_2.fastq.gz" \
            "$DATADIR/${s}_1.fastq.gz" "$DATADIR/${s}_2.fastq.gz" \
            > "$TRIM_DIR/${s}.log" 2>&1

        wait_for_file "$TRIM_DIR/${s}_clean_1.fastq.gz" "trim $s R1"
        wait_for_file "$TRIM_DIR/${s}_clean_2.fastq.gz" "trim $s R2"
    else
        cutadapt -q 20 -m 25 \
            -o "$TRIM_DIR/${s}_clean.fastq.gz" \
            "$DATADIR/${s}.fastq.gz" \
            > "$TRIM_DIR/${s}.log" 2>&1

        wait_for_file "$TRIM_DIR/${s}_clean.fastq.gz" "trim $s SE"
    fi
done
