#!/usr/bin/env bash
set -euo pipefail

DATADIR="${DATADIR:-$1}"
TRIM_DIR="${TRIM_DIR:-$2}"
mkdir -p "$TRIM_DIR"

echo "[1/6] Trimming FASTQ files in $DATADIR → $TRIM_DIR"

files=( "$DATADIR"/*.fastq.gz )
[[ ${#files[@]} -gt 0 ]] || { echo "[ERROR] No FASTQ files in $DATADIR"; exit 1; }

declare -A samples
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

for sample in "${!samples[@]}"; do
    layout="${samples[$sample]}"
    if [[ "$layout" == "PE" ]]; then
        r1="$DATADIR/${sample}_1.fastq.gz"
        r2="$DATADIR/${sample}_2.fastq.gz"
        clean_r1="$TRIM_DIR/${sample}_clean_1.fastq.gz"
        clean_r2="$TRIM_DIR/${sample}_clean_2.fastq.gz"
        echo "[PE] $sample"
        "$CUTADAPT" -q 20 -m 25 -o "$clean_r1" -p "$clean_r2" "$r1" "$r2" > "$TRIM_DIR/${sample}_cutadapt.log" 2>&1
    else
        r1="$DATADIR/${sample}.fastq.gz"
        clean_se="$TRIM_DIR/${sample}_clean.fastq.gz"
        echo "[SE] $sample"
        "$CUTADAPT" -q 20 -m 25 -o "$clean_se" "$r1" > "$TRIM_DIR/${sample}_cutadapt.log" 2>&1
    fi
done

echo "[OK] Trimming completed"
