#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ==================================== USAGE DESCRIPTION ==================================== #
usage() {
    cat <<EOF
Usage: ./concatenate_fastq.sh -m manifest -f
    -m, --manifest path to the TSV file
    -c, --fast_cat Concatenate .gz without recompressing
    -f, --force Overwrite outputs if they already exist
    -h, --help Shows help
The TSV file must contain the columns (tab): sample_name layout R1 R2 concat_dir out_prefix
EOF
}
# =========================================================================================== #

# ==================================== PARSING ARGUMENTS ==================================== #
MANIFEST=""
FORCE=false
FAST_CAT=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--manifest) MANIFEST="$2"; shift 2 ;;
        -f|--force)    FORCE=true; shift ;;
        -c|--fast_cat) FAST_CAT=true; shift ;;
        -h|--help)     usage; exit 0 ;;
        --) shift; break ;;
        -*) echo "Opción inválida: $1" >&2; usage; exit 2 ;;
        *)  break ;;
    esac
done
# =========================================================================================== #

# ==================================== VERIFICATIONS ==================================== #
[[ -z "${MANIFEST}" ]] && { echo "[ERROR] -m|--manifest is absent" >&2; usage; exit 2; }
[[ -r "${MANIFEST}" ]] || { echo "[ERROR] ${MANIFEST} couldn't be read" >&2; exit 2; }
# ======================================================================================= #

# ==================================== FUNCTIONS ==================================== #
COMPRESS_CMD=(gzip -c)

check_fastq() { # Validates FASTQ files (existence, size and structure)
    local f="$1"
    # -s evaluates
    [[ -s "${f}" ]] && gzip -t "${f}" >/dev/null 2>&1 || return 1
    zcat "${f}" | awk 'NR%4==1 && substr($0,1,1)!="@"{b=1}
                        NR%4==3 && substr($0,1,1)!="+"{c=1}
                        END{exit b||c}' || return 1
    return 0
}

count_reads() {
  zcat "$1" | awk 'END{print (NR/4)}'
}

concat_list() {
  # $1 = file list (string with \n), $2 = output
  local list="$1" out="$2"
  mkdir -p "$(dirname -- "${out}")"
  local tmp="${out}.tmp.$$"

  # Building array from list
  mapfile -t files <<< "${list}"

  # If the array is empty, then it is skipped
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "[INFO] No hay archivos para ${out}. Saltando."
    return 0
  fi

  # Validate inputs (existence, gzip, FASTQ expected structure)
  for f in "${files[@]}"; do
    [[ -n "${f}" ]] || continue
    check_fastq "${f}" || { echo "[ERROR] Invalid FASTQ: ${f}" >&2; rm -f "${tmp}"; return 1; }
  done

  # Output already existing
  if [[ -f "${out}" && "${FORCE}" != "true" ]]; then
    echo "[INFO] ${out} file already exists. Use --force to overwrite. Skipping."
    return 0
  fi

  if [[ "${FAST_CAT}" == "true" ]]; then
    echo "[INFO] Concatenating (fast_cat, without recompression) -> ${out}"
    cat "${files[@]}" > "${tmp}"
  else
    echo "[INFO] Concatenating (gzip recompression) -> ${out}"
    zcat "${files[@]}" | "${COMPRESS_CMD[@]}" > "${tmp}"
  fi

  # Validating output file (complete gzip)
  gzip -t "${tmp}" >/dev/null 2>&1 || { echo "[ERROR] Corrupt output: ${tmp}" >&2; rm -f "${tmp}"; return 1; }

  mv -f "${tmp}" "${out}"
  return 0
}
# =================================================================================== #

# ==================================== READING MANIFEST ==================================== #
declare -A R1_MAP
declare -A R2_MAP

# Read TSV (ignores empty lines and comments with #)
while IFS=$'\t' read -r sample_name layout R1 R2 concat_dir out_prefix; do
  # Skip empty lines
  [[ -z "${sample_name:-}" ]] && continue
  # Skip comments
  [[ "${sample_name:0:1}" == "#" ]] && continue
  # Skip common header
  if [[ "${sample_name}" == "sample_name" && "${layout}" == "layout" ]]; then
    continue
  fi

  key="${concat_dir}|${out_prefix}|${layout}|${sample_name}"

  # Acumulate lists per group (separated by a newline)
  R1_MAP["$key"]+="${R1}"$'\n'
  if [[ "${layout}" == "PAIRED" && -n "${R2:-}" && "${R2}" != "-" ]]; then
    R2_MAP["$key"]+="${R2}"$'\n'
  fi
done < "${MANIFEST}"
# ========================================================================================== #

# =================== PROCESSING EACH GROUP =================== #
for key in "${!R1_MAP[@]}"; do
  IFS='|' read -r concat_dir out_prefix layout sample_name <<< "${key}"
  # Standarizing final path
  concat_dir="${concat_dir%/}"
  out_r1="${concat_dir}/${out_prefix}_1.fastq.gz"
  out_r2="${concat_dir}/${out_prefix}_2.fastq.gz"

  # Concat R1
  concat_list "${R1_MAP[$key]}" "${out_r1}"

  if [[ "${layout}" == "PAIRED" ]]; then
    # Concat R2
    concat_list "${R2_MAP[$key]:-}" "${out_r2}"

    # Parity R1 vs R2
    if [[ -f "${out_r1}" && -f "${out_r2}" ]]; then
      n1=$(count_reads "${out_r1}")
      n2=$(count_reads "${out_r2}")
      if [[ "${n1}" -ne "${n2}" ]]; then
        echo "[WARNING] Not equal lectures for ${sample_name} (${out_prefix}): R1=${n1}, R2=${n2}" >&2
      else
        echo "[INFO] Correct lectures parity for ${sample_name} (${out_prefix}): ${n1} lecturas."
      fi
    fi
  fi

  echo "[DONE] ${sample_name} → ${out_prefix} en ${concat_dir}"
done

# ============================================================ #