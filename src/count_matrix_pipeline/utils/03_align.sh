#!/usr/bin/env bash
set -euo pipefail
# NOTE: intended to be sourced

echo "=== [3/6] Running BWA-MEM alignment ==="

BWA_INDEX="$bwa_index_prefix"

for sample in "${!samples[@]}"; do
    layout="${samples[$sample]}"
    sam_out="$BWA_DIR/${sample}.sam"
    log_out="$BWA_DIR/${sample}_bwa.log"

    echo "--- Processing $sample ($layout) ---"

    if [[ "$layout" == "PE" ]]; then
        r1="$TRIM_DIR/${sample}_clean_1.fastq.gz"
        r2="$TRIM_DIR/${sample}_clean_2.fastq.gz"

        if [[ -f "$r1" && -f "$r2" ]]; then
            echo "[PE] Aligning $sample"
            echo "Input: $r1, $r2"
            echo "Output: $sam_out"

            bwa mem -t 4 "$BWA_INDEX" "$r1" "$r2" > "$sam_out" 2> "$log_out"
            exit_code=$?

            # Wait for SAM to be ready
            wait_for_file "$sam_out"
        else
            echo "ERROR: Missing trimmed files for $sample"
            exit 1
        fi
    else
        r1="$TRIM_DIR/${sample}_clean.fastq.gz"

        if [[ -f "$r1" ]]; then
            echo "[SE] Aligning $sample"
            echo "Input: $r1"
            echo "Output: $sam_out"

            bwa mem -t 8 "$BWA_INDEX" "$r1" > "$sam_out" 2> "$log_out"
            exit_code=$?

            wait_for_file "$sam_out"
        else
            echo "ERROR: Missing trimmed file for $sample"
            exit 1
        fi
    fi

    # Verify file
    if [[ $exit_code -eq 0 && -s "$sam_out" ]]; then
        lines=$(wc -l < "$sam_out")
        echo "SUCCESS: $sample - $lines lines in SAM"
    else
        echo "FAILED: $sample - exit code: $exit_code"
        echo "BWA log output:"
        cat "$log_out" || true
        echo "SAM file exists: $(if [[ -f "$sam_out" ]]; then echo "YES ($(wc -l < "$sam_out") lines)"; else echo "NO"; fi)"
        exit 1
    fi
    echo
done

echo "[OK] All BWA alignments completed successfully"
