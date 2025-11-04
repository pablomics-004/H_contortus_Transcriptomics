#!/usr/bin/env bash
set -euo pipefail
BWA="${BWA:-bwa}"
REF_GENOME="${REF_GENOME:-$1}"
BWA_PREFIX="${BWA_PREFIX:-$2}"
BWA_DIR="${BWA_DIR:-$3}"
FORCE_INDEX="${FORCE_INDEX:-0}"

echo "[2/6] Creating BWA index"

BWA_INDEX_DIR="$BWA_DIR/BWA_INDEX"
mkdir -p "$BWA_INDEX_DIR"
bwa_index_prefix="$BWA_INDEX_DIR/$BWA_PREFIX"

[[ -s "$REF_GENOME" ]] || { echo "[ERROR] Reference genome $REF_GENOME not found"; exit 1; }

if [[ ! -f "${bwa_index_prefix}.bwt" || $FORCE_INDEX -eq 1 ]]; then
    echo "[INFO] Building BWA index..."
    "$BWA" index -p "$bwa_index_prefix" "$REF_GENOME" > "${bwa_index_prefix}_index.log" 2>&1

    for ext in amb ann bwt pac sa; do
        idx_file="${bwa_index_prefix}.${ext}"
        timeout=300
        waited=0
        while [[ ! -f "$idx_file" ]]; do
            sleep 5
            waited=$((waited+5))
            if [[ $waited -ge $timeout ]]; then
                echo "[ERROR] Timeout waiting for BWA index file $idx_file"
                exit 1
            fi
        done
    done
else
    echo "[OK] BWA index already exists"
fi

echo "[OK] All BWA index files ready"
