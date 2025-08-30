#!/usr/bin/env bash
set -euo pipefail

# Directory with the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DOC="$SCRIPT_DIR/../doc"
DATA="$SCRIPT_DIR/../data"
TMP="$SCRIPT_DIR/../tmp"
RESULTS="$SCRIPT_DIR/../results"

mkdir -p "$DOC" "$DATA" "$TMP" "$RESULTS"

shopt -s nullglob

links=(
    "$DOC"
    "$DATA"
    "$RESULTS"
)

outfile="$TMP/tmp_file.txt"
printf "Temporal file to let the visualization" >> "$outfile"

for d in "${links[@]}"; do
    ln -s "$outfile" "$d"
done

shopt -u nullglob