#!/usr/bin/env bash
set -euo pipefail

run="${1:?Required run name}"
run_dir="${2:?Required run directory}"
paired_end="${3:-1}"
thr_fq="${4:-4}"
thr_pz="${5:-6}"


[[ ! -d "$run_dir" ]] && echo "The directory ${run_dir} doesn't exist" && exit 1
[[ $thr_fq -gt 8 ]] && echo "[WARNING] It's not recommendable to use more than 8 threads within the server"
[[ $thr_pz -gt 8 ]] && echo "[WARNING] It's not recommendable to use more than 8 threads within the server"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "[ERROR] Required command '$1' not found in PATH" >&2
    exit 127
  }
}

check_fastq_format() {
    local fq_file="$1"
    awk '
        NR % 4 == 1 && $0 !~ /^@/ {
            bad = 1; print "[ERROR] Line " NR ": expected @ at beginning";
            exit 1 
        }
        NR % 4 == 3 && $0 !~ /^\+/ {
            bad = 1; print "[ERROR] Line " NR ": expected + at beginning";
            exit 1
        }
        END {
            if (NR % 4 != 0) {
                print "[ERROR] File length not multiple of 4 (" NR " lines)";
                exit 1
            }
        }
    ' "$fq_file" || {
        echo "[FAIL] Invalid FASTQ format: ${fq_file}" >&2
        exit 1
    }
}

# require_cmd prefetch it's not necessary because fasterq-dump retrieves directly from NCBI
require_cmd fasterq-dump
require_cmd pigz


echo "[INFO] Running fasterq-dump on $run..."
fasterq-dump "$run" --threads $thr_fq --split-files --outdir "$run_dir"

fq1="${run_dir}/${run}_1.fastq"
fq2="${run_dir}/${run}_2.fastq"

echo "[BASH] Checking FASTQ format for ${fq1}..."
check_fastq_format "$fq1"

if (( paired_end )); then
    echo "[BASH] Checking FASTQ format for ${fq2}..."
    check_fastq_format "$fq2"
fi

echo "[BASH] Compressing ${run}_1.fastq..."
pigz -p $thr_pz -1 "${run_dir}/${run}_1.fastq"

if (( paired_end )); then
    echo "[BASH] Compressing ${run}_2.fastq..."
    pigz -p $thr_pz -1 "${run_dir}/${run}_2.fastq"
fi

echo "[BASH] Completed successfully for $run"