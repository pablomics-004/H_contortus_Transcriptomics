#!/usr/bin/env bash
set -euo pipefail
# NOTE: intended to be sourced

echo "=== [2/6] Creating BWA index ==="

BWA_INDEX_DIR="$BWA_DIR/BWA_INDEX"
mkdir -p "$BWA_INDEX_DIR"
bwa_index_prefix="$BWA_INDEX_DIR/$BWA_PREFIX"

[[ -s "$REF_GENOME" ]] || { echo "[ERROR] Reference genome $REF_GENOME not found"; exit 1; }

if [[ ! -f "${bwa_index_prefix}.bwt" ]]; then
    echo "[INFO] Building BWA index..."
    bwa index -p "$bwa_index_prefix" "$REF_GENOME" > "${bwa_index_prefix}_index.log" 2>&1

    # Wait for index files (use wait_for_file on at least one; then check all)
    echo "Waiting for BWA index files..."
    wait_for_file "${bwa_index_prefix}.bwt"
else
    echo "[OK] BWA index already exists."
fi

# Verify index files
index_files=("${bwa_index_prefix}.amb" "${bwa_index_prefix}.ann" "${bwa_index_prefix}.bwt" "${bwa_index_prefix}.pac" "${bwa_index_prefix}.sa")
for index_file in "${index_files[@]}"; do
    if [[ ! -f "$index_file" ]]; then
        echo "ERROR: BWA index file missing: $(basename "$index_file")"
        exit 1
    fi
done
echo "[OK] All BWA index files ready"

