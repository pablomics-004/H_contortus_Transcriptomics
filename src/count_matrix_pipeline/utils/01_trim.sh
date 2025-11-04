#!/usr/bin/env bash
set -euo pipefail
# NOTE: this script is intended to be sourced by launch_pipe.sh
# so variables like DATADIR, TRIM_DIR are expected to be set.

echo "=== [1/6] Starting Cutadapt trimming ==="

files=( "$DATADIR"/*.fastq.gz )
[[ ${#files[@]} -gt 0 ]] || { echo "[ERROR] No FASTQ files found in $DATADIR"; exit 1; }

# Identificar samples únicos (PE/SE)
declare -gA samples
for f in "${files[@]}"; do
    base=$(basename "$f")
    if [[ "$base" == *_1.fastq.gz ]]; then
        sample="${base%_1.fastq.gz}"
        samples["$sample"]="PE"
    elif [[ "$base" == *_2.fastq.gz ]]; then
        continue
    else
        sample="${base%.fastq.gz}"
        samples["$sample"]="SE"
    fi
done

echo "Found ${#samples[@]} samples: ${!samples[*]}"

# Ejecutar cutadapt directamente (sin qsub), uno por uno y esperar salidas
for sample in "${!samples[@]}"; do
    layout="${samples[$sample]}"

    if [[ "$layout" == "PE" ]]; then
        r1="$DATADIR/${sample}_1.fastq.gz"
        r2="$DATADIR/${sample}_2.fastq.gz"
        clean_r1="$TRIM_DIR/${sample}_clean_1.fastq.gz"
        clean_r2="$TRIM_DIR/${sample}_clean_2.fastq.gz"

        if [[ -f "$r1" && -f "$r2" ]]; then
            echo "[PE] Processing $sample"
            cutadapt -q 20 -m 25 -o "$clean_r1" -p "$clean_r2" "$r1" "$r2" \
                > "$TRIM_DIR/${sample}_cutadapt.log" 2>&1

            # Wait until both files exist and are non-empty
            wait_for_file "$clean_r1"
            wait_for_file "$clean_r2"

            # Double-check exit status from last command: cutadapt exit code is in $?
            if [[ ${PIPESTATUS[0]:-0} -eq 0 || $? -eq 0 ]]; then
                echo "$sample: PE trimming successful"
            else
                echo "$sample: PE trimming failed"
                cat "$TRIM_DIR/${sample}_cutadapt.log" || true
                exit 1
            fi
        else
            echo "[ERROR] Missing files for $sample: $r1 or $r2"
            exit 1
        fi
    else
        r1="$DATADIR/${sample}.fastq.gz"
        clean_se="$TRIM_DIR/${sample}_clean.fastq.gz"

        if [[ -f "$r1" ]]; then
            echo "[SE] Processing $sample"
            cutadapt -q 20 -m 25 -o "$clean_se" "$r1" \
                > "$TRIM_DIR/${sample}_cutadapt.log" 2>&1

            wait_for_file "$clean_se"

            if [[ ${PIPESTATUS[0]:-0} -eq 0 || $? -eq 0 ]]; then
                echo "$sample: SE trimming successful"
            else
                echo "$sample: SE trimming failed"
                cat "$TRIM_DIR/${sample}_cutadapt.log" || true
                exit 1
            fi
        else
            echo "[ERROR] Missing file for $sample: $r1"
            exit 1
        fi
    fi
done

# Verify trimmed files (redundant given wait_for_file, but kept as final check)
for sample in "${!samples[@]}"; do
    if [[ "${samples[$sample]}" == "PE" ]]; then
        if [[ ! -f "$TRIM_DIR/${sample}_clean_1.fastq.gz" || ! -f "$TRIM_DIR/${sample}_clean_2.fastq.gz" ]]; then
            echo "[ERROR] Trimmed files missing for $sample"
            exit 1
        fi
    else
        if [[ ! -f "$TRIM_DIR/${sample}_clean.fastq.gz" ]]; then
            echo "[ERROR] Trimmed file missing for $sample"
            exit 1
        fi
    fi
done
echo "[OK] All trimmed files verified"
