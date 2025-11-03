#!/usr/bin/env bash
set -euo pipefail

thr_pz="${1:-6}"
outfile="${2:?Required outfile name}"
shift 2

if (( $# == 0 )); then
    echo "[ERROR] No input files provided to concatenate" >&2
    exit 1
fi

echo "[BASH] Concatenating ${#@} files into $outfile"
zcat "$@" | pigz -p "$thr_pz" -1 > "$outfile"