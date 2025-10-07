#!/usr/bin/env bash
set -euo pipefail

# Usage: ./prefet_concaten.sh sra_by_biosample.csv /path/sra_dir /path/concat_dir
csv="${1:?CSV file required}"
sra_dir="${2:?sra_dir required}"
concat_dir="${3:?concat_dir required}"

# ----------------------------- Helpers & Guards ------------------------------

# Require a command to exist or exit with a clear message
require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "[error] Required command '$1' not found in PATH" >&2
    exit 127
  }
}

require_cmd prefetch
require_cmd fastq-dump
require_cmd gzip
# zcat is generally available via gzip package; check as well:
require_cmd zcat || true

mkdir -p "$sra_dir" "$concat_dir"

# Concatenate .gz → .gz with pigz if available (faster), else gzip
concat_gz() {
  local out="$1"; shift
  if command -v pigz >/dev/null 2>&1; then
    echo "[*] Compressing with pigz → $out"
    pigz -dc "$@" | pigz > "$out"
  else
    echo "[*] Compressing with gzip → $out"
    zcat "$@" | gzip > "$out"
  fi
}

# Decompress a single .gz to stdout using pigz if available
decompress() {
  if command -v pigz >/devnull 2>&1; then
    pigz -dc "$1"
  else
    zcat "$1"
  fi
}

# Validate a .fastq.gz file:
#  - exists and is non-empty
#  - gzip container passes integrity test
#  - FASTQ block structure holds:
#      * line 1 of each record starts with '@'
#      * line 3 of each record starts with '+'
#      * len(line 4) == len(line 2)
check_fastq_gz() {
  local f="$1"
  [[ -s "$f" ]] || { echo "[warn] Empty or missing FASTQ: $f"; return 1; }
  gzip -t "$f" >/dev/null 2>&1 || { echo "[warn] Corrupted gzip: $f"; return 1; }
  decompress "$f" | awk '
    NR%4==1 { if (substr($0,1,1)!="@") bad=1; next }
    NR%4==2 { sl=length($0); next }
    NR%4==3 { if (substr($0,1,1)!="+") bad=1; next }
    NR%4==0 { if (length($0)!=sl) bad=1; next }
    END { exit bad }
  ' || { echo "[warn] FASTQ structure invalid: $f"; return 1; }
  return 0
}

# ----------------------------- CSV → TSV (awk) -------------------------------

# CSV parsing that respects quoted fields (handles "SRR1,SRR2")
awk_cmd='
BEGIN { FPAT = "([^,])|(\"([^\"]|\"\")\")" }
NR>1 {
  biosample = $1; srr = $2; layout = $3
  gsub(/^"|"$/, "", biosample)
  gsub(/^"|"$/, "", srr)
  gsub(/^"|"$/, "", layout)
  print biosample "\t" srr "\t" layout
}'

# --------------------------------- Main Loop ---------------------------------

while IFS=$'\t' read -r biosample srr_list layout; do
  [[ -z "$biosample" ]] && continue

  safe_sample=$(printf "%s" "$biosample" | tr -c '[:alnum:]' '_')
  srr_list="${srr_list//,/ }"   # "SRR1,SRR2" -> "SRR1 SRR2"

  r1_files=()
  r2_files=()

  echo "[*] Processing BioSample: $biosample (layout: ${layout:-UNKNOWN})"

  for run in $srr_list; do
    [[ -z "$run" ]] && continue

    run_dir="$sra_dir/$run"
    mkdir -p "$run_dir"

    r1_gz="$run_dir/${run}_1.fastq.gz"
    r2_gz="$run_dir/${run}_2.fastq.gz"

    echo "[*] Prefetch $run"
    prefetch -O "$sra_dir" "$run"

    echo "[*] fastq-dump $run"
    if [[ "${layout^^}" == "PAIRED" ]]; then
      # For paired-end, emit _1/_2 and gzip directly
      fastq-dump "$run" -O "$run_dir" --split-files --gzip
    else
      # For single-end, fastq-dump emits ${run}.fastq.gz; normalize to *_1.fastq.gz
      fastq_gz="$run_dir/${run}.fastq.gz"
      if [[ ! -s "$fastq_gz" && ! -s "$r1_gz" ]]; then
        fastq-dump "$run" -O "$run_dir" --gzip
      fi
      [[ -s "$fastq_gz" && ! -e "$r1_gz" ]] && mv "$fastq_gz" "$r1_gz"
    fi

    # ---- Per-run FASTQ validation ----
    if check_fastq_gz "$r1_gz"; then
      echo "[ok] R1 validated: $r1_gz"
      r1_files+=("$r1_gz")
    else
      echo "[warn] R1 invalid for $run → skipping this run"
      continue
    fi

    if [[ "${layout^^}" == "PAIRED" ]]; then
      if [[ -s "$r2_gz" ]]; then
        if check_fastq_gz "$r2_gz"; then
          echo "[ok] R2 validated: $r2_gz"
          r2_files+=("$r2_gz")
        else
          echo "[warn] R2 invalid for $run → skipping this run"
          # remove the R1 we added for this run to keep pairs consistent
          unset 'r1_files[${#r1_files[@]}-1]'
          continue
        fi
      else
        echo "[warn] Missing R2 for $run → skipping this run"
        # remove the R1 we added for this run to keep pairs consistent
        unset 'r1_files[${#r1_files[@]}-1]'
        continue
      fi
    fi
  done

  # ------------------------ Per-sample concatenation -------------------------

  if [[ ${#r1_files[@]} -gt 0 ]]; then
    out1="$concat_dir/${safe_sample}_1.fastq.gz"
    echo "[*] Concatenate R1 → $out1"
    concat_gz "$out1" "${r1_files[@]}"
    echo "[*] Validate concatenated R1"
    if ! check_fastq_gz "$out1"; then
      echo "[error] Invalid concatenated R1: $out1" >&2
      exit 1
    fi
  else
    echo "[warn] No valid R1 files for sample $biosample"
  end

  if [[ "${layout^^}" == "PAIRED" && ${#r2_files[@]} -gt 0 ]]; then
    out2="$concat_dir/${safe_sample}_2.fastq.gz"
    echo "[*] Concatenate R2 → $out2"
    concat_gz "$out2" "${r2_files[@]}"
    echo "[*] Validate concatenated R2"
    if ! check_fastq_gz "$out2"; then
      echo "[error] Invalid concatenated R2: $out2" >&2
      exit 1
    fi
  elif [[ "${layout^^}" == "PAIRED" ]]; then
    echo "[warn] No valid R2 files for paired sample $biosample"
  fi

  echo "[ok] Sample completed: $biosample"
  echo
done < <(awk -F, "$awk_cmd" "$csv")