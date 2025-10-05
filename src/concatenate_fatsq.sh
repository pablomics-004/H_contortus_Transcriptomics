#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
    cat <<EOF
Usage: ./concatenate_fastq.sh -m manifest -f
    -m, --manifest path to the TSV file
    -f, --force Overwrite outputs if they already exist
    -h, --help Shows help
The TSV file must contain the columns (tab): sample_name layout R1 R2 concat_dir out_prefix
EOF
}

MANIFEST=""
FORCE=false

while [[ $# -gt 0 ]]; do
    case "$$1" in
        -m|--manifest) MANIFEST="$2"; shift 2 ;;
        -f|--force)    FORCE=true; shift ;;
        -h|--help)     usage; exit 0 ;;
        --) shift; break ;;
        -*) echo "Opción inválida: $1" >&2; usage; exit 2 ;;
        *)  break ;;
    esac
done

[[ -z "${MANIFEST}" ]] && { echo "[ERROR] -m|--manifest is absent" >&2; usage; exit 2; }
[[ -r "${MANIFEST}" ]] || { echo "[ERROR] ${MANIFEST} couldn't be read" >&2; exit 2; }

COMPRESS_CMD=("gzip -c")

check_fastq() { # Verifies if it is valid
    local f="$1"
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
    local list="$1" out="$2"
    mkdir -p "$(dirname -- "${out}")"
    local tmp="${out}.tmp.$$"

    # String to an array
    mapfile -t files <<< "${list}"

    [[ ${#files[@]} -eq 0 ]] && { echo "[ERROR] Empty list for ${out}" >&2; return 1; }


    for f in "${files[@]}"; do
        [[ -n "${f}" ]] || continue
        check_fastq "${f}" || { echo "[ERROR] FASTQ inválido: ${f}" >&2; rm -f "${tmp}"; return 1; }
    done

    if [[ -f "${out}" && "${FORCE}" != "true" ]]; then
        echo "[INFO] File ${out} already exists. Use --force to overwrite. Skipping."
        return 0
    fi

    echo "[INFO] Concatenating -> ${out}"
    zcat "${files[@]}" | "${COMPRESS_CMD[@]}" > "${tmp}"
    gzip -t "${tmp}" >/dev/null 2>&1
    mv -f "${tmp}" "${out}"
    return 0
}

# Pairwise arrays, works like a python dictionary
declare -A R1_MAP
declare -A R2_MAP

while IFS=$'\t' read -r sample_name layout R1 R2 concat_dir out_prefix; do
  [[ -z "${sample_name:-}" ]] && continue
  [[ "${sample_name:0:1}" == "#" ]] && continue

  key="${concat_dir}|${out_prefix}|${layout}|${sample_name}"

  R1_MAP["$key"]+="${R1}"$'\n'
  if [[ "${layout}" == "PAIRED" && -n "${R2:-}" ]]; then
    R2_MAP["$key"]+="${R2}"$'\n'
  fi
done < "${MANIFEST}"

for key in "${!R1_MAP[@]}"; do
  IFS='|' read -r concat_dir out_prefix layout sample_name <<< "${key}"
  out_r1="${concat_dir%/}/${out_prefix}_1.fastq.gz"
  out_r2="${concat_dir%/}/${out_prefix}_2.fastq.gz"

  concat_list "${R1_MAP[$key]}" "${out_r1}"

  if [[ "${layout}" == "PAIRED" ]]; then
    concat_list "${R2_MAP[$key]:-}" "${out_r2}"

    if [[ -f "${out_r1}" && -f "${out_r2}" ]]; then
      n1=$(count_reads "${out_r1}")
      n2=$(count_reads "${out_r2}")
      if [[ "${n1}" -ne "${n2}" ]]; then
        echo "[WARN] Desbalance en lecturas para ${sample_name} (${out_prefix}): R1=${n1}, R2=${n2}" >&2
      else
        echo "[INFO] Paridad OK para ${sample_name} (${out_prefix}): ${n1} lecturas."
      fi
    fi
  fi

  echo "[DONE] ${sample_name} → ${out_prefix} en ${concat_dir}"
done