"""
BioProject Entrez Pipeline
--------------------------
A modular pipeline for querying, downloading, and organizing data from NCBI BioProject,
including BioSample metadata, SRA runs, and reference genomes.

Author: Ashley Yael Montiel Vargas & Pablo Salazar Mendez
Version: 1.0
Date: 2025-10-07

This script integrates multiple NCBI Entrez databases using Biopython, automating:
 - Retrieval of BioProject metadata
 - Identification of linked BioSamples and SRA accessions
 - Download and concatenation of FASTQ files from SRA
 - Retrieval of reference genome sequences and annotations

Dependencies:
 - Biopython (Entrez module)
 - pandas
 - subprocess, urllib, xml.etree.ElementTree
 - NCBI SRA Toolkit (prefetch, fastq-dump)

Example usage:
    python3 bioproject_pipeline.py -p PRJNA877658 -e user@email.com -i -b -s --concat -r
"""

from pathlib import Path
from Bio import Entrez
from pandas import DataFrame
import logging
import subprocess
import os
from collections import defaultdict
import urllib.request
from datetime import datetime
import xml.etree.ElementTree as ET
import argparse

# ==================== LOGGING CONFIGURATION ==================== #
log_filename = f"bioproject_pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(log_filename), logging.StreamHandler()]
)
# =============================================================== #

# ==================== UTILITIES ==================== #
def file_exists(filepath: str) -> bool:
    """
    Check whether a file exists and is non-empty.
    """
    return os.path.exists(filepath) and os.path.getsize(filepath) > 0

def check_format(filepath: str, expected_extension: str) -> bool:
    """
    Validate that a file matches the expected format based on its extension.
    Logs a warning if mismatched.
    """

    if not filepath.endswith(expected_extension):
        logging.warning(f"Format mismatch: {filepath} (expected {expected_extension})")
        return False
    return True

def concat_gz(files: list[str], output_file: str):
    """
    Concatenate multiple compressed FASTQ files (.fastq.gz) into one output file.

    Uses `zcat` and `gzip` via subprocess for efficiency with compressed data.
    Skips concatenation if output file already exists.
    """
    if not files:
        return
    if file_exists(output_file):
        logging.info(f"{output_file} already exists. Skipping concatenation.")
        return
    logging.info(f"Concatenating {len(files)} files to {output_file}")
    try:
        with open(output_file, "wb") as fout:
            p1 = subprocess.Popen(["zcat"] + files, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=fout)
            p1.stdout.close()
            p2.communicate()
        if check_format(output_file, ".fastq.gz"):
            logging.info(f"Concatenation completed: {output_file}")
        else:
            logging.warning(f"Format mismatch after concatenation: {output_file}")
    except Exception as e:
        logging.error(f"Error during concatenation: {e}")

# ==================== ENTIREZ CONFIG ==================== #
def configure_entrez(email: str, api_key: str | None = None):
    """
    Configure Entrez API credentials.

    Parameters:
        email (str): User email (required by NCBI policy).
        api_key (str, optional): NCBI API key for higher request limits.
    """
    if not email:
        raise ValueError("Email is required for NCBI Entrez")
    Entrez.email = email
    Entrez.api_key = api_key if api_key else None
    logging.info(f"Entrez configured with email: {email}, API key: {'Provided' if api_key else 'None'}")
    print(f"Entrez configured with email: {email}")


# ==================== BIOPROJECT INFO ==================== #
def get_bioproject_info(bioproject_id: str) -> tuple[dict, str] | None:
    """
    Retrieve general metadata for a given BioProject accession.
    Returns (dict, str): BioProject summary and internal UID.
    """
    try:
        logging.info(f"Fetching BioProject: {bioproject_id}")
        with Entrez.esearch(db="bioproject", term=bioproject_id) as handle:
            record = Entrez.read(handle)
        if not record["IdList"]:
            logging.warning(f"No BioProject found: {bioproject_id}")
            return None
        uid = record["IdList"][0]
        with Entrez.esummary(db="bioproject", id=uid) as summary:
            results = Entrez.read(summary)
        result = results["DocumentSummarySet"]["DocumentSummary"][0]
        msg = (
            f"\n=== BioProject Information ===\n"
            f"Title: {result.get('Project_Title', 'N/A')}\n"
            f"Accession: {result.get('Project_Acc', 'N/A')}\n"
            f"Description: {result.get('Project_Description', 'N/A')}\n"
            f"Data Type: {result.get('Project_Data_Type', 'N/A')}\n"
            f"Submitter: {result.get('Submitter_Organization', 'N/A')}\n"
            f"Registration Date: {result.get('Registration_Date', 'N/A')}\n"
            f"URL: https://www.ncbi.nlm.nih.gov/bioproject/{bioproject_id}\n"
        )
        logging.info(msg)
        return result, uid
    except Exception as e:
        logging.error(f"Error fetching BioProject info: {e}")
        return None

# ==================== DATABASE LINKS ==================== #
def get_associated_databases(uid: str):
     """
    Retrieve a list of databases associated with a given BioProject UID.

    Returns:
        list[str]: Sorted list of linked database names.
    """
    try:
        with Entrez.elink(dbfrom="bioproject", id=uid) as handle:
            records = Entrez.read(handle)
        linked_dbs = sorted({
            db['DbTo']
            for linkset in records if 'LinkSetDb' in linkset
            for db in linkset['LinkSetDb']
        })
        if linked_dbs:
            logging.info("=== Associated Databases ===\n" + "\n".join(f"- {db}" for db in linked_dbs))
        else:
            logging.info(f"No associated databases found for {uid}")
        return linked_dbs
    except Exception as e:
        logging.error(f"Error retrieving associated databases: {e}")
        return []

# ==================== BIOSAMPLES ==================== #
def fetch_biosamples(bioproject_uid: str, output: str = None):
    """
    Retrieve all BioSamples linked to a given BioProject UID.

    Parameters:
        bioproject_uid (str): Internal UID of the BioProject.
        output (str, optional): CSV file path to save the results.

    Returns:
        pandas.DataFrame: BioSample metadata.
    """
    try:
        with Entrez.elink(dbfrom='bioproject', id=bioproject_uid, db='biosample') as handle:
            records = Entrez.read(handle)
        biosample_ids = [
            link['Id']
            for linkset in records
            for db in linkset.get('LinkSetDb', [])
            for link in db.get('Link', [])
            if 'Id' in link
        ]
        if not biosample_ids:
            logging.warning('No BioSamples linked to this BioProject')
            return None
        results = []
        for sample_id in biosample_ids:
            try:
                with Entrez.esummary(db='biosample', id=sample_id) as handle:
                    summary = Entrez.read(handle)["DocumentSummarySet"]["DocumentSummary"][0]
                identifiers = summary.get("Identifiers", "")
                sra = identifiers.split("SRA:")[-1].strip() if "SRA:" in identifiers else "N/A"
                results.append({
                    'BioSample': sample_id,
                    'Accession': summary.get('Accession', 'N/A'),
                    'Organism': summary.get('Organism', 'N/A'),
                    'Title': summary.get('Title', 'N/A'),
                    'SubmissionDate': summary.get('Date', 'N/A'),
                    'SRA': sra,
                    'Attributes': summary.get('Attributes', 'N/A')
                })
            except Exception as e:
                logging.warning(f"Error fetching BioSample {sample_id}: {e}")
        df = DataFrame(results).drop_duplicates(subset=["Accession"])
        if output:
            df.to_csv(output, index=False)
            logging.info(f"BioSamples saved to {output}")
        return df
    except Exception as e:
        logging.error(f"Error fetching BioSamples: {e}")
        return None

# ==================== SRA INFO ==================== #
def get_sra_by_bioproject(bioproject_uid: str, output_csv: str):
    """
    Retrieve all SRA accessions (SRR) and their corresponding BioSamples from a given BioProject.

    Parameters:
        bioproject_uid (str): BioProject UID.
        output_csv (str): Path to save the resulting CSV file.

    Returns:
        pandas.DataFrame: Mapping of BioSamples to SRR accessions and library layouts.
    """
    sra_by_sample = defaultdict(list)
    layout_by_sample = {}
    with Entrez.elink(dbfrom="bioproject", db="sra", id=bioproject_uid) as handle:
        linkset = Entrez.read(handle)
    sra_ids = {
        link['Id']
        for db in linkset[0].get('LinkSetDb', [])
        for link in db.get('Link', [])
    }
    logging.info(f"Found {len(sra_ids)} SRA entries linked to BioProject {bioproject_uid}")

    for sra_id in sra_ids:
        try:
            with Entrez.efetch(db="sra", id=sra_id, rettype="xml") as handle:
                xml_data = handle.read()
            root = ET.fromstring(xml_data)
            for exp_pkg in root.findall(".//EXPERIMENT_PACKAGE"):
                sample_node = exp_pkg.find(".//SAMPLE")
                biosample = sample_node.get("accession", sample_node.get("alias", "UnknownSample")) if sample_node is not None else "UnknownSample"
                layout_node = exp_pkg.find(".//LIBRARY_LAYOUT")
                if layout_node is not None:
                    layout = "PAIRED" if layout_node.find("PAIRED") is not None else "SINGLE"
                else:
                    layout = "PAIRED"
                layout_by_sample[biosample] = layout
                for run in exp_pkg.findall(".//RUN"):
                    srr = run.get("accession")
                    if srr:
                        sra_by_sample[biosample].append(srr)
                        sra_by_sample[biosample] = list(dict.fromkeys(sra_by_sample[biosample]))
        except Exception as e:
            logging.warning(f"Error fetching SRA {sra_id}: {e}")

    df = DataFrame([
        {"BioSample": sample, "SRR_List": ",".join(srrs), "LibraryLayout": layout_by_sample[sample]}
        for sample, srrs in sra_by_sample.items()
    ])
    df.to_csv(output_csv, index=False)
    logging.info(f"SRA info saved to {output_csv}")
    return df

# ==================== REFERENCE GENOME ==================== #
def get_reference_genome_info(organism: str, outdir="reference_genome", max_results=5):
    """
    Retrieve the latest reference genome for a given organism from NCBI Assembly.

    Parameters:
        organism (str): Scientific name of the organism (e.g., "Homo sapiens").
        outdir (str, optional): Directory to save the downloaded files. Default is "reference_genome".
        max_results (int, optional): Maximum number of assemblies to retrieve. Default is 5.

    Returns:
        None
    """
    os.makedirs(outdir, exist_ok=True)
    logging.info(f"Searching reference genome for: {organism}")
    try:
        with Entrez.esearch(db="assembly", term=f"{organism}[organism] AND latest[filter] AND reference[filter]", retmax=max_results) as handle:
            record = Entrez.read(handle)
        if not record["IdList"]:
            logging.warning(f"No reference genome found for {organism}")
            return
        for asm_id in record["IdList"]:
            with Entrez.esummary(db="assembly", id=asm_id) as summary:
                doc = Entrez.read(summary)
                info_list = doc.get("DocumentSummarySet", {}).get("DocumentSummary", [])
                if not info_list:
                    continue
                info = info_list[0]
            ftp_path = info.get("FtpPath_RefSeq") or info.get("FtpPath_GenBank")
            if not ftp_path:
                continue
            asm_name = ftp_path.split("/")[-1]
            files_to_download = [
                (f"{ftp_path}/{asm_name}_genomic.fna.gz", f"{asm_name}_genomic.fna.gz"),
                (f"{ftp_path}/{asm_name}_genomic.gff.gz", f"{asm_name}_genomic.gff.gz")
            ]
            for url, filename in files_to_download:
                filepath = os.path.join(outdir, filename)
                if file_exists(filepath):
                    logging.info(f"{filename} exists, skipping.")
                    continue
                try:
                    logging.info(f"Downloading {filename} ...")
                    urllib.request.urlretrieve(url, filepath)
                    logging.info(f"{filename} downloaded successfully")
                except Exception as e:
                    logging.error(f"Error downloading {filename}: {e}")
    except Exception as e:
        logging.error(f"Error fetching reference genome: {e}")

# ==================== SRA DOWNLOAD & CONCAT ==================== #

def download_and_concat_sra_grouped(sra_df: DataFrame, sra_dir: str, concat_dir: str):
    """
    Download all SRA runs grouped by BioSample and concatenate FASTQ files per sample.

    Parameters:
        sra_df (pandas.DataFrame): DataFrame containing BioSample names, SRR accessions, and library layouts.
        sra_dir (str): Directory to store individual downloaded SRA runs.
        concat_dir (str): Directory to store concatenated FASTQ files per sample.

    Returns:
        None
    """
    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(concat_dir, exist_ok=True)
    grouped = defaultdict(list)
    layout_map = {}
    for _, row in sra_df.iterrows():
        sample_name = row['BioSample']
        runs = row['SRR_List'].split(",")
        grouped[sample_name].extend(runs)
        layout_map[sample_name] = row.get('LibraryLayout', 'PAIRED')

    for sample_name, runs in grouped.items():
        r1_files, r2_files = [], []
        for run in runs:
            run_dir = os.path.join(sra_dir, run)
            os.makedirs(run_dir, exist_ok=True)
            r1_path = os.path.join(run_dir, f"{run}_1.fastq.gz")
            r2_path = os.path.join(run_dir, f"{run}_2.fastq.gz") if layout_map[sample_name] == "PAIRED" else None
            if not file_exists(r1_path) or (r2_path and not file_exists(r2_path)):
                logging.info(f"Prefetching {run} ...")
                try:
                    subprocess.run(["prefetch", "--output-directory", sra_dir, run], check=True)
                    logging.info(f"Converting {run} to FASTQ.gz ...")
                    subprocess.run(["fastq-dump", run, "-O", run_dir, "--split-files", "--gzip"], check=True)
                except subprocess.CalledProcessError as e:
                    logging.error(f"Error downloading {run}: {e}")
                    continue
            if file_exists(r1_path):
                r1_files.append(r1_path)
            if r2_path and file_exists(r2_path):
                r2_files.append(r2_path)

        safe_sample_name = "".join(c if c.isalnum() else "_" for c in sample_name)
        if r1_files:
            concat_gz(r1_files, os.path.join(concat_dir, f"{safe_sample_name}_1.fastq.gz"))
        if layout_map[sample_name] == "PAIRED" and r2_files:
            concat_gz(r2_files, os.path.join(concat_dir, f"{safe_sample_name}_2.fastq.gz"))
        logging.info(f"Sample {sample_name} processed successfully.")

# ==================== MAIN ==================== #
def main():
    """
    Main function to run the BioProject pipeline.

    This function parses command-line arguments to query and download data from NCBI Entrez, including BioProject info, 
    associated databases, BioSamples, SRA runs, and reference genomes. It also handles downloading and concatenating 
    FASTQ files per BioSample.

    Command-line arguments:
        -p / --project: BioProject ID (required)
        -e / --email: Email for NCBI Entrez (required)
        -a / --api_key: Optional API key for NCBI Entrez
        -o / --outdir: Output directory (default: "bioproject_data")
        -i / --info: Show BioProject information
        -d / --databases: Show associated databases
        -b / --biosamples: Retrieve BioSamples information
        -s / --samples: Generate SRA samples table
        --download-all: Download all SRA FASTQ runs
        --concat: Concatenate multiple FASTQ files per sample
        -r / --reference: Download reference genome (optional)
        --max-ref-results: Maximum reference genomes to retrieve (default: 5)

    Returns:
        None
    """

    parser = argparse.ArgumentParser(description="Pipeline to query and download BioProject data using Entrez.")
    parser.add_argument("-p", "--project", required=True, help="BioProject ID (e.g., PRJNA877658)")
    parser.add_argument("-e", "--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("-a", "--api_key", help="API key for NCBI Entrez (optional)")
    parser.add_argument("-o", "--outdir", default="bioproject_data", help="Output directory")
    parser.add_argument("-i", "--info", action="store_true", help="Show BioProject info")
    parser.add_argument("-d", "--databases", action="store_true", help="Show associated databases")
    parser.add_argument("-b", "--biosamples", action="store_true", help="Retrieve BioSamples info")
    parser.add_argument("-s", "--samples", action="store_true", help="Get SRA samples table")
    parser.add_argument("--download-all", action="store_true", help="Download all SRA FASTQ runs")
    parser.add_argument("--concat", action="store_true", help="Concatenate multiple FASTQs per sample")
    parser.add_argument("-r", "--reference", nargs='?', const="auto", metavar="ORGANISM", help="Download reference genome")
    parser.add_argument("--max-ref-results", type=int, default=5, help="Max reference genomes to show")
    parser.add_argument("--use-bash",action="store_true",help="Use the external Bash script (prefet_concaten.sh) for SRA download & concatenation instead of Python")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sra_dir = os.path.join(args.outdir, "sra")
    concat_dir = os.path.join(args.outdir, "fastq_merged")
    biosample_dir = os.path.join(args.outdir, "biosamples")
    ref_dir = os.path.join(args.outdir, "reference_genome")
    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(concat_dir, exist_ok=True)
    os.makedirs(biosample_dir, exist_ok=True)
    os.makedirs(ref_dir, exist_ok=True)

    configure_entrez(args.email, args.api_key)
    logging.info(f"Pipeline started for {args.project}")

    project_info = None
    if args.info:
        project_info, uid = get_bioproject_info(args.project)
    if args.databases:
        get_associated_databases(uid if project_info else args.project)
    df_biosamples = None
    if args.biosamples:
        df_biosamples = fetch_biosamples(uid, output=os.path.join(biosample_dir, "biosamples.csv"))
    df_sra = None
    if args.samples or args.download_all or args.concat:
        df_sra = get_sra_by_bioproject(uid, output_csv=os.path.join(sra_dir, "sra_by_biosample.csv"))
    if args.download_all or args.concat:
        if df_sra is not None and not df_sra.empty:
            if args.use_bash:
                # Prepare paths
                csv_path = os.path.join(sra_dir, "sra_by_biosample.csv")
                df_sra.to_csv(csv_path, index=False)

                bash_script = "./prefet_concaten.sh"  # Adjust if path differs
                logging.info(f"Running Bash SRA pipeline: {bash_script}")
                try:
                    subprocess.run(
                        [bash_script, csv_path, sra_dir, concat_dir],
                        check=True
                    )
                    logging.info("Bash SRA pipeline completed successfully.")
                except subprocess.CalledProcessError as e:
                    logging.error(f"Bash SRA pipeline failed: {e}")
            else:
                # Python fallback
                download_and_concat_sra_grouped(df_sra, sra_dir, concat_dir)
    if args.reference:
        organism_to_search = args.reference
        if organism_to_search == "auto" and df_biosamples is not None:
            organism_to_search = df_biosamples['Organism'].drop_duplicates().tolist()[0]
        if organism_to_search:
            get_reference_genome_info(organism_to_search, outdir=ref_dir, max_results=args.max_ref_results)

    logging.info("Pipeline execution completed successfully.")
    print(f"\nPipeline completed! Log file: {log_filename}")

if __name__ == "__main__":
    main()
