#!/usr/bin/env python3
# ============================================================
# bioproject_pipeline.py
# Pipeline para consultar, descargar y procesar datos de BioProject (NCBI)
# ============================================================

import os
import sys
import subprocess
import argparse
import pandas as pd
import urllib.request
import logging
from Bio import Entrez

# ==================== CONFIGURACIÓN ==================== #
logging.basicConfig(filename="bioproject_pipeline.log",
                    level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")
# ======================================================== #

# ==================== UTILIDADES ==================== #
def file_exists(filepath: str) -> bool:
    """Verifica si un archivo existe y no está vacío."""
    return os.path.exists(filepath) and os.path.getsize(filepath) > 0

def check_format(filepath: str, expected_extension: str) -> bool:
    """Verifica que el archivo tenga la extensión esperada."""
    if not filepath.endswith(expected_extension):
        logging.warning(f"Format mismatch: {filepath} (expected {expected_extension})")
        return False
    return True
# ======================================================== #

# ==================== SRA HANDLING ==================== #
def download_sra_run(srr_id: str, outdir: str):
    """Descarga un SRR y lo convierte a FASTQ.gz"""
    os.makedirs(outdir, exist_ok=True)
    logging.info(f"Downloading {srr_id} to {outdir}")

    subprocess.run(["prefetch", "--output-directory", outdir, srr_id], check=True)
    srr_path = os.path.join(outdir, srr_id)

    subprocess.run(["fastq-dump", "--split-files", "--gzip", "-O", srr_path, srr_id], check=True)
    logging.info(f"Downloaded and extracted {srr_id}")

def concatenate_fastqs(sample_name: str, files_r1: list, files_r2: list = None, outdir: str = "."):
    """Concatena archivos FASTQ.gz (paired o single)."""
    os.makedirs(outdir, exist_ok=True)
    out_r1 = os.path.join(outdir, f"{sample_name}_1.fastq.gz")

    print(f"Concatenating {len(files_r1)} R1 files → {out_r1}")
    with open(out_r1, "wb") as fout:
        p1 = subprocess.Popen(["zcat"] + files_r1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=fout)
        p1.stdout.close()
        p2.communicate()

    if files_r2:
        out_r2 = os.path.join(outdir, f"{sample_name}_2.fastq.gz")
        print(f"Concatenating {len(files_r2)} R2 files → {out_r2}")
        with open(out_r2, "wb") as fout:
            p1 = subprocess.Popen(["zcat"] + files_r2, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=fout)
            p1.stdout.close()
            p2.communicate()

    logging.info(f"Concatenation completed for {sample_name}")
# ======================================================== #

# ==================== REFERENCE GENOME ==================== #
def get_reference_genome_info(organism: str, outdir="reference_genome", max_results=5):
    """Descarga el genoma de referencia y anotaciones (FASTA + GFF)."""
    os.makedirs(outdir, exist_ok=True)
    print(f"Searching reference genome for: {organism}")
    try:
        with Entrez.esearch(db="assembly",
                            term=f"{organism}[organism] AND latest[filter] AND reference[filter]",
                            retmax=max_results) as handle:
            record = Entrez.read(handle)
        if not record["IdList"]:
            print(f"No reference genome found for {organism}")
            return
        for i, asm_id in enumerate(record["IdList"]):
            with Entrez.esummary(db="assembly", id=asm_id) as summary:
                doc = Entrez.read(summary)
                info = doc["DocumentSummarySet"]["DocumentSummary"][0]
            ftp_path = info.get("FtpPath_RefSeq") or info.get("FtpPath_GenBank")
            if not ftp_path:
                print(f"No FTP path found for assembly {asm_id}")
                continue
            asm_name = ftp_path.split("/")[-1]
            print(f"\n--- Reference Assembly {i+1} ---")
            print(f"Name: {info.get('AssemblyName', 'N/A')}")
            print(f"Accession: {info.get('AssemblyAccession', 'N/A')}")
            print(f"Species: {info.get('SpeciesName', 'N/A')}")
            print(f"FTP URL: {ftp_path}")

            fasta_url = f"{ftp_path}/{asm_name}_genomic.fna.gz"
            gff_url = f"{ftp_path}/{asm_name}_genomic.gff.gz"

            for url, filename in [(fasta_url, f"{asm_name}_genomic.fna.gz"),
                                  (gff_url, f"{asm_name}_genomic.gff.gz")]:
                filepath = os.path.join(outdir, filename)
                if file_exists(filepath):
                    print(f"{filename} already exists. Skipping.")
                    continue
                try:
                    print(f"Downloading {filename} ...")
                    urllib.request.urlretrieve(url, filepath)
                    if file_exists(filepath) and check_format(filepath, ".gz"):
                        print(f"Successfully downloaded {filename}")
                    else:
                        print(f"Failed to download {filename}")
                except Exception as e:
                    logging.error(f"Error downloading {filename}: {e}")
        print(f"Reference genome download completed for {organism}")
    except Exception as e:
        logging.error(f"Error fetching reference genome: {e}")
# ======================================================== #

# ==================== MAIN PIPELINE ==================== #
def main():
    parser = argparse.ArgumentParser(description="Pipeline to query and download BioProject data using Entrez.")
    parser.add_argument("-p", "--project", required=True, help="BioProject ID (e.g., PRJNA877658)")
    parser.add_argument("-e", "--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("-a", "--api_key", help="API key for NCBI Entrez (optional)")
    parser.add_argument("-o", "--outdir", default="bioproject_data", help="Output directory for all downloads")
    parser.add_argument("--download-all", action="store_true", help="Download all SRA FASTQ runs")
    parser.add_argument("--concat", action="store_true", help="Concatenate multiple FASTQs per sample")
    parser.add_argument("-r", "--reference", nargs='?', const="auto", metavar="ORGANISM",
                        help="Download reference genome. Without argument uses BioProject organism, or specify an organism")
    parser.add_argument("--max-ref-results", type=int, default=5,
                        help="Maximum number of reference results to show (default: 5)")
    args = parser.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    os.makedirs(args.outdir, exist_ok=True)

    # Ejemplo: lista de SRRs simulada (en práctica se obtiene con Entrez)
    srr_list = ["SRR031714", "SRR031715"]

    if args.download_all:
        for srr in srr_list:
            download_sra_run(srr, os.path.join(args.outdir, "sra_runs"))

    if args.concat:
        sample_name = "sampleGSM461177"
        base_dir = os.path.join(args.outdir, "sra_runs")
        r1_files, r2_files = [], []
        for srr in srr_list:
            r1 = os.path.join(base_dir, srr, f"{srr}_1.fastq.gz")
            r2 = os.path.join(base_dir, srr, f"{srr}_2.fastq.gz")
            if file_exists(r1):
                r1_files.append(r1)
            if file_exists(r2):
                r2_files.append(r2)
        if r1_files:
            concatenate_fastqs(sample_name, r1_files, r2_files if r2_files else None, outdir=args.outdir)

    if args.reference:
        organism = args.reference if args.reference != "auto" else "Haemonchus contortus"
        get_reference_genome_info(organism, outdir=os.path.join(args.outdir, "reference_genome"),
                                  max_results=args.max_ref_results)

if __name__ == "__main__":
    main()
