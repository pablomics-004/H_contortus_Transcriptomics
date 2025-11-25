#!/usr/bin/env bash

###############################################################################
# Script: process_trinotate.sh
# Description:
#     Processes a gene list and a Trinotate annotation matrix to extract:
#         1. KO terms
#         2. GO terms
#         3. Filter annotations by %ID > 55 or helminth lineage
#
# Input:
#     1) gene_list.txt            (1 column, gene IDs)
#     2) trinotate_report.xls     (exactly 14 tab-separated columns)
#
# Output:
#     genes_annotation.tsv
#     filtered_genes.tsv
#     genes_kos_gos.tsv
#     DE_genes_KO.tsv
#     DE_genes_GO.tsv
#
# Usage:
#     ./process_trinotate.sh gene_list.txt trinotate.xls
###############################################################################

set -euo pipefail

###########################
# Validate arguments
###########################
if [[ $# -ne 2 ]]; then
    echo "Error: This script requires exactly two input files:"
    echo "  1) gene_list.txt"
    echo "  2) Trinotate report (14-column TSV/XLS)"
    exit 1
fi

GENE_LIST="$1"
TRINOTATE="$2"

if [[ ! -f "$GENE_LIST" ]]; then
    echo "Error: Gene list file not found: $GENE_LIST"
    exit 1
fi
if [[ ! -f "$TRINOTATE" ]]; then
    echo "Error: Trinotate file not found: $TRINOTATE"
    exit 1
fi

#############################################
# Validate gene list: must have 1 column
#############################################
cols_gene_list=$(awk '{print NF; exit}' "$GENE_LIST")
if [[ $cols_gene_list -ne 1 ]]; then
    echo "Error: gene_list.txt must contain exactly 1 column."
    echo "Detected columns: $cols_gene_list"
    exit 1
fi

#############################################
# Validate Trinotate file: must have 14 columns
#############################################
cols_trinotate=$(awk -F"\t" '{print NF; exit}' "$TRINOTATE")
if [[ $cols_trinotate -ne 14 ]]; then
    echo "Error: Trinotate file must contain exactly 14 tab-separated columns."
    echo "Detected columns: $cols_trinotate"
    exit 1
fi

echo "Input validation completed successfully."
echo "---------------------------------------------"


###############################################################################
# STEP 1 — Keep only genes present in the gene list
###############################################################################
echo "Step 1 — Extracting annotation only for listed genes..."

awk '
    NR==FNR { genes[$1]; next }
    $2 in genes
' "$GENE_LIST" "$TRINOTATE" > genes_annotation.tsv

echo "Generated: genes_annotation.tsv"
echo


###############################################################################
# STEP 2 — Filter by %ID > 55 or helminth taxonomy
###############################################################################
echo "Step 2 — Filtering by percent identity and helminth lineage..."

awk -F"\t" '
function getID(s) {
    match(s, /[0-9.]+%ID/, a)
    gsub("%ID","",a[0])
    return a[0] + 0
}
function isHelminth(s) {
    return (s ~ /Nematoda/ || s ~ /Platyhelminthes/ || s ~ /Acanthocephala/)
}

NR==1 { print; next }

{
    id = getID($3)
    if (id > 55 || isHelminth($3))
        print
}
' genes_annotation.tsv > filtered_genes.tsv

echo "Generated: filtered_genes.tsv"
echo


###############################################################################
# STEP 3 — Keep rows with at least one KO or GO term
###############################################################################
echo "Step 3 — Retaining rows with KO or GO annotations..."

awk -F"\t" '
NR==1 { print; next }

{
    kegg = $12
    gob  = $13
    gop  = $14

    if ((kegg != "." && kegg != "") ||
        (gob != "." && gob != "") ||
        (gop != "." && gop != ""))
        print
}
' filtered_genes.tsv > genes_kos_gos.tsv

echo "Generated: genes_kos_gos.tsv"
echo


###############################################################################
# STEP 4 — Extract KO terms
###############################################################################
echo "Step 4 — Extracting KO terms..."

awk -F'\t' '
$12 != "." {
    s = $12
    while(match(s, /KO:K[0-9]+/)) {
        print $2 "\t" substr(s, RSTART+3, RLENGTH-3)
        s = substr(s, RSTART+RLENGTH)
    }
}
' genes_kos_gos.tsv > DE_genes_KO.tsv

echo "Generated: DE_genes_KO.tsv"
echo


###############################################################################
# STEP 5 — Extract GO terms
###############################################################################
echo "Step 5 — Extracting GO terms..."

awk -F"\t" '
{
    for(i=13;i<=14;i++) {
        if($i != ".") {
            s = $i
            while(match(s, /GO:[0-9]+/)) {
                print $2 "\t" substr(s, RSTART, RLENGTH)
                s = substr(s, RSTART+RLENGTH)
            }
        }
    }
}
' genes_kos_gos.tsv > DE_genes_GO.tsv

echo "Generated: DE_genes_GO.tsv"
echo
echo "Pipeline completed successfully."
