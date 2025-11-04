#!/usr/bin/env bash
log "=== [2/6] BWA Index ==="

[[ -s "$REF" ]] || die "Reference fasta empty"

prefix="$BWA_DIR/BWA_INDEX/$BWA_PREFIX"

if [[ ! -f "${prefix}.bwt" ]]; then
    bwa index -p "$prefix" "$REF" > "${prefix}.log" 2>&1
    wait_for_file "${prefix}.bwt" "BWA index"
else
    log "Index already exists"
fi
