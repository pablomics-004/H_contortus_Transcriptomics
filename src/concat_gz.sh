#!/usr/bin/env bash
set -euo pipefail

outfile="${1:?Required outfile name}"; shift

if (( $# == 0 )); then
    echo "[ERROR] No input files provided to concatenate" >&2
    exit 1
fi

echo "[INFO] Concatenating ${#@} files into $outfile"

zcat "$@" | pigz -p 6 -1 > "$outfile"