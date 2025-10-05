#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
    cat <<EOF
Usage: ./src/pf_fastqc.sh -s <sra_paths> [-p|--paired] SRR21518936 SRR21518939 ...
    -s, --sra-dir   Output directory (subdirectories creation for each SRR)
    -c, --concat-dir    Directory for concatenation per sample
    -S, --sample-name   Sample name/label (Title)
    -o, --out-prefix    Output prefix per sample (by default = sample_name)
    -m, --manifest      Generated TSV path (default: <sra_dir>/manifest_runs.tsv)
    -h, --help  Shows help to the user

Example:
    ./src/pf_fastqc.sh -s /data/sra -p SRR21518936 SRR21518939 ...
EOF
}

# Initializing variables
PAIRED=false
SRA_DIR=""
SAMPLE_NAME=""
OUT_PREFIX=""
MANIFEST=""

# Parsing the arguments of the script
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--sra-dir) SRA_DIR="$2"; shift 2 ;;
    -c|--concat-dir)  CONCAT_DIR="$2"; shift 2 ;;
    -S|--sample-name) SAMPLE_NAME="$2"; shift 2 ;;
    -o|--out-prefix)  OUT_PREFIX="$2"; shift 2 ;;
    -m|--manifest)    MANIFEST="$2"; shift 2 ;;
    -p|--paired)  PAIRED=true; shift ;;
    -h|--help)    usage; exit 0 ;;
    --) shift; break ;;
    -*) echo "Invalid option: $1" >&2; usage; exit 2 ;;
    *)  break ;;
  esac
done

[[ -z "${SRA_DIR}"      ]] && { echo "[ERROR] -s|--sra-dir is absent" >&2; usage; exit 2; }
[[ -z "${CONCAT_DIR}"   ]] && { echo "[ERROR] -c|--concat-dir  is absent" >&2; usage; exit 2; }
[[ -z "${SAMPLE_NAME}"  ]] && { echo "[ERROR] -S|--sample-name  is absent" >&2; usage; exit 2; }
[[ $# -lt 1             ]] && { echo "[ERROR] provide at least one SRR" >&2; usage; exit 2; }

OUT_PREFIX="${OUT_PREFIX:-${SAMPLE_NAME}}"
MANIFEST="${MANIFEST:-${SRA_DIR%/}/manifest_runs.tsv}"
LAYOUT=$($PAIRED && echo "PAIRED" || echo "SINGLE")

mkdir -p "${SRA_DIR}" "${CONCAT_DIR}"

# If the manifest doesn't exists
if [[ ! -s "${MANIFEST}" ]]; then
  printf "# sample_name\tlayout\tR1\tR2\tconcat_dir\tout_prefix\n" > "${MANIFEST}"
fi


have_fastq() {
    local run_dir="$1" srr="$2" paired="$3"
    local r1="${run_dir}/${srr}_1.fastq.gz"
    local r2="${run_dir}/${srr}_2.fastq.gz"

    # Checking out the existence and integrity of the first file
    [[ -s "${r1}" ]] && gzip -t "${r1}" >/dev/null 2>&1 || return 1
    zcat "${r1}" | awk 'NR%4==1 && substr($0,1,1)!="@"{b=1}
                        NR%4==3 && substr($0,1,1)!="+"{c=1}
                        END{exit b||c}' || return 1

    if [[ "${paired}" == "true" ]]; then
        [[ -s "${r2}" ]] && gzip -t "${r2}" >/dev/null 2>&1 || return 1
        zcat "${r2}" | awk 'NR%4==1 && substr($0,1,1)!="@"{b=1}
                        NR%4==3 && substr($0,1,1)!="+"{c=1}
                        END{exit b||c}' || return 1
    fi
    return 0
}

for SRR in "$@"; do
    run_dir="${SRA_DIR}/${SRR}"
    mkdir -p "${run_dir}"

    echo "[INFO] SRR=${SRR}"
    if have_fastq "${run_dir}" "${SRR}" "${PAIRED}"; then
        echo "  - FASTQ files present and intact. Skipping."
        continue
    fi

    echo "  - Prefetch ${SRR} ..."
    prefetch --output-directory "${SRA_DIR}" "${SRR}"

    echo "  - Converting to FASTQ.gz (${SRR}) ..."
    fastq-dump "${SRR}" -O "${run_dir}" --split-files --gzip

    # Validating output files
    have_fastq "${run_dir}" "${SRR}" "${PAIRED}" || { echo "  - [WARN] FASTQ inválidos: ${SRR}" >&2; exit 1; }
  else
    echo "  - FASTQ presentes y válidos. Saltando."
  fi
    # Final paths
    r1="${run_dir}/${SRR}_1.fastq.gz"
    r2=""
    if [[ "${LAYOUT}" == "PAIRED" ]]; then
    r2="${run_dir}/${SRR}_2.fastq.gz"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${SAMPLE_NAME}" "${LAYOUT}" "${r1}" "${r2}" "${CONCAT_DIR}" "${OUT_PREFIX}" >> "${MANIFEST}"
done

echo "[DONE] Download and conversion completed."