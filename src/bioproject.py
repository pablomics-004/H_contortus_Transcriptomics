"""
bioproject_pipeline.py
----------------------------------
Comprehensive BioProject pipeline using Entrez.

Features:
- Retrieve BioProject, BioSample, and SRA information.
- Confirm valid project access.
- Download and validate sequence files (FASTQ/FASTA/GFF).
- Concatenate paired-end or single-end files automatically without intermediate files.
- Retrieve and download reference genomes from Assembly (RefSeq).
- Logging and modular structure with argparse.

Author: Ashley Yael Montiel Vargas & Pablo Salazar Mendez
Date: Oct-2025
"""

from Bio import Entrez
import pandas as pd
import os
import subprocess
import argparse as ap
import logging
import time
import urllib.request
from collections import defaultdict
from datetime import datetime

# ==================== LOGGING CONFIGURATION ==================== #
log_filename = f"bioproject_pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler()
    ]
)
# =============================================================== #

# ==================== ENTIREZ CONFIG ==================== #
def configure_entrez(email: str, api_key: str | None = None):
    if not email:
        raise ValueError("Email is required for NCBI Entrez")
    Entrez.email = email
    Entrez.api_key = api_key if api_key else None
    logging.info(f"Entrez configured with email: {email}, API key: {'Provided' if api_key else 'None'}")
    print(f"Entrez configured with email: {email}")
# ========================================================== #

# ==================== UTILITY ==================== #
def file_exists(filepath: str) -> bool:
    return os.path.exists(filepath) and os.path.getsize(filepath) > 0

def check_format(filepath: str, expected_extension: str) -> bool:
    if not filepath.endswith(expected_extension):
        logging.warning(f"Format mismatch: {filepath} (expected {expected_extension})")
        return False
    return True
# =================================================== #

# ==================== BIOPROJECT ==================== #
def get_bioproject_info(bioproject_id: str) -> dict | None:
    try:
        with Entrez.esearch(db="bioproject", term=f"{bioproject_id}[project_acc]") as handle:
            record = Entrez.read(handle)
        if not record["IdList"]:
            raise ValueError(f"Couldn't find BioProject {bioproject_id}")
        uid = record["IdList"][0]
        with Entrez.esummary(db="bioproject", id=uid) as summary:
            result = Entrez.read(summary)[0]
        msg = (f"\n=== BioProject Information ===\n"
               f"Title: {result.get('Project_Title', 'N/A')}\n"
               f"Accession: {result.get('Project_Acc', 'N/A')}\n"
               f"Description: {result.get('Project_Description', 'N/A')}\n"
               f"Organism: {result.get('Organism_Name', 'N/A')}\n"
               f"Submission Date: {result.get('Submission_Date', 'N/A')}\n"
               f"Last Update: {result.get('Last_Update', 'N/A')}\n"
               f"URL: https://www.ncbi.nlm.nih.gov/bioproject/{bioproject_id}\n")
        print(msg)
        logging.info(msg)
        return result
    except Exception as e:
        logging.error(f"Error fetching BioProject info: {e}")
        print(f"Error fetching BioProject info: {e}")
        return None
# =================================================== #

# ==================== DATABASE LINKS ==================== #
def get_associated_databases(bioproject_id: str) -> list[str]:
    try:
        with Entrez.elink(dbfrom="bioproject", id=bioproject_id) as handle:
            records = Entrez.read(handle)
        linked_dbs = []
        for linkset in records:
            if "LinkSetDb" in linkset:
                for db in linkset["LinkSetDb"]:
                    linked_dbs.append(db["DbTo"])
        linked_dbs = sorted(set(linked_dbs))
        msg = "=== Associated Databases ===\n" + "\n".join(f"- {db}" for db in linked_dbs)
        print(msg)
        logging.info(msg)
        return linked_dbs
    except Exception as e:
        logging.error(f"Error retrieving associated databases: {e}")
        print(f"Error retrieving associated databases: {e}")
        return []
# =================================================== #

# ================= FETCH BIOSAMPLES ================= #
def fetch_biosamples(bid: str) -> dict | None:
    try:
        with Entrez.esummary(db='biosample', id=bid) as handle:
            summary = Entrez.read(handle)[0]
        return {
            'BioSample_ID' : summary.get('Id', dft := 'N/A'),
            'Accession' : summary.get('Accession', dft),
            'Organism' : summary.get('Organism', dft),
            'Title' : summary.get('Title', dft),
            'Attributes' : summary.get('Attributes', dft),
            'SubmissionDate' : summary.get('SubmissionDate', dft)
        }
    except Exception as inner_e:
        logging.warning(f"Error retrieving BioSample {bid}: {inner_e}")
        return
    finally:
        time.sleep(0.34)
# =================================================== #

# ==================== BIOSAMPLES ==================== #
def get_biosample_info(bioproject_id: str, output: str) -> pd.DataFrame | None:
    try:
        with Entrez.elink(dbfrom="bioproject", id=bioproject_id, db="biosample") as handle:
            records = Entrez.read(handle)
        
        # biosample_ids = []
        # for linkset in records:
        #     for db in linkset.get("LinkSetDb", []):
        #         for link in db.get("Link", []):
        #             biosample_ids.append(link["Id"])
        
        biosample_ids = [
            link['Id']
            for linkset in records
            for db in linkset.get('LinkSetDb', [])
            for link in db.get('Link', []) if 'Id' in link
        ]

        if not biosample_ids:
            print("No BioSamples linked to this BioProject.")
            return
                
        # biosamples = []
        # for bid in biosample_ids:
        #     try:
        #         with Entrez.esummary(db="biosample", id=bid) as handle:
        #             summary = Entrez.read(handle)[0]
        #         biosamples.append({
        #             "BioSample_ID": summary.get("Id", "N/A"),
        #             "Accession": summary.get("Accession", "N/A"),
        #             "Organism": summary.get("Organism", "N/A"),
        #             "Title": summary.get("Title", "N/A"),
        #             "Attributes": summary.get("Attributes", "N/A"),
        #             "SubmissionDate": summary.get("SubmissionDate", "N/A")
        #         })
        #         time.sleep(0.3)
        #     except Exception as inner_e:
        #         logging.warning(f"Error retrieving BioSample {bid}: {inner_e}")
        
        df = pd.DataFrame(filter(None, map(fetch_biosamples, biosample_ids)))
        df.to_csv(output, index=False)
        print(f"Retrieved {len(df)} BioSamples. Saved to {output}")
        logging.info(f"Retrieved {len(df)} BioSamples.")
        return df
    
    except Exception as e:
        logging.error(f"Error fetching BioSamples: {e}")
        print(f"Error fetching BioSamples: {e}")
        return
# =================================================== #

# ================== FETCH SRA ================== #
def fetch_sra(sra_id: str) -> dict | None:
    try:
        with Entrez.esummary(db="sra", id=sra_id) as handle:
            summary = Entrez.read(handle)[0]
            return {
                'RunAccession' : summary.get('Runs', dft:='N/A').split(',')[0],
                'Title' : summary.get('Title', dft),
                'Organism' : summary.get('Organism', dft),
                'Instrument' : summary.get('Instrument', dft),
                'LibraryLayout' : summary.get('LibraryLayout', dft),
                'Bases' : summary.get('Bases', dft),
                'Spots' : summary.get('Spots', dft)
            }
    except Exception as inner_e:
        logging.warning(f"Skipping SRA ID {sra_id}: {inner_e}")
        return
    finally:
        time.sleep(0.34)
# =============================================== #

# ==================== SRA ==================== #
def get_sra_samples(bioproject_id: str, output: str) -> pd.DataFrame | None:
    try:
        with Entrez.esearch(db="sra", term=f"{bioproject_id}[bioproject]", retmax=500) as handle:
            record = Entrez.read(handle)
        ids = record["IdList"]
        print(f"Found {len(ids)} SRA entries linked to {bioproject_id}")
        logging.info(f"Found {len(ids)} SRA entries linked to {bioproject_id}")
        if not ids:
            print("Warning: No SRA entries found.")
            return
        
        # sample_data = []
        # for sra_id in ids:
        #     try:
        #         with Entrez.esummary(db="sra", id=sra_id) as handle:
        #             summary = Entrez.read(handle)[0]
        #         sample_data.append({
        #             "RunAccession": summary.get("Runs", "N/A").split(",")[0],
        #             "Title": summary.get("Title", "N/A"),
        #             "Organism": summary.get("Organism", "N/A"),
        #             "Instrument": summary.get("Instrument", "N/A"),
        #             "LibraryLayout": summary.get("LibraryLayout", "N/A"),
        #             "Bases": summary.get("Bases", "N/A"),
        #             "Spots": summary.get("Spots", "N/A")
        #         })
        #         time.sleep(0.3)
        #     except Exception as inner_e:
        #         logging.warning(f"Skipping SRA ID {sra_id}: {inner_e}")

        df = pd.DataFrame(filter(None, map(fetch_sra, ids)))
        df.to_csv(output, index=False)
        print(f"Saved SRA sample info to {output}")
        logging.info(f"Saved {len(df)} SRA entries to {output}")
        return df
    except Exception as e:
        logging.error(f"Error fetching SRA samples: {e}")
        print(f"Error fetching SRA samples: {e}")
        return
# =================================================== #

# ==================== SRA DOWNLOAD & CONCAT ==================== #
def download_and_concat_sra_grouped(sra_df: pd.DataFrame, sra_dir: str, concat_dir: str):

    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(concat_dir, exist_ok=True)

    grouped = defaultdict(list)
    layout_map = {}

    for _, row in sra_df.iterrows():
        grouped[row['Title']].append(row['RunAccession'])
        layout_map[row['Title']] = row['LibraryLayout']

    for sample_name, runs in grouped.items():
        r1_files, r2_files = [], []

        for run in runs:
            run_dir = os.path.join(sra_dir, run)
            os.makedirs(run_dir, exist_ok=True)
            r1_path = os.path.join(run_dir, f"{run}_1.fastq.gz")
            r2_path = os.path.join(run_dir, f"{run}_2.fastq.gz") if layout_map[sample_name] == "PAIRED" else None

            if not file_exists(r1_path) or (r2_path and not file_exists(r2_path)):
                print(f"Prefetching {run} ...")
                subprocess.run(["prefetch", "--output-directory", sra_dir, run], check=True)
                print(f"Converting {run} to FASTQ.gz ...")
                subprocess.run(["fastq-dump", run, "-O", run_dir, "--split-files", "--gzip"], check=True)
            else:
                print(f"{run} already downloaded. Skipping.")

            if file_exists(r1_path) and check_format(r1_path, ".fastq.gz"):
                r1_files.append(r1_path)
            if r2_path and file_exists(r2_path) and check_format(r2_path, ".fastq.gz"):
                r2_files.append(r2_path)

        def concat_gz(files, output_file):
            if file_exists(output_file):
                print(f"{output_file} already exists. Skipping concatenation.")
                return
            with open(output_file, "wb") as fout:
                p1 = subprocess.Popen(["zcat"] + files, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=fout)
                p1.stdout.close()
                p2.communicate()
            if check_format(output_file, ".fastq.gz"):
                print(f"Concatenation completed: {output_file}")
            else:
                print(f"Warning: format mismatch after concatenation: {output_file}")

        if r1_files:
            out_r1 = os.path.join(concat_dir, f"{sample_name}_1.fastq.gz")
            concat_gz(r1_files, out_r1)
        if layout_map[sample_name] == "PAIRED" and r2_files:
            out_r2 = os.path.join(concat_dir, f"{sample_name}_2.fastq.gz")
            concat_gz(r2_files, out_r2)

        print(f"Sample {sample_name} processed successfully.\n")
        
# =================================================== #

# ==================== REFERENCE GENOME ==================== #
def get_reference_genome_info(organism: str, outdir="reference_genome", max_results=5) -> None:
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

            files_to_download = [
                (f"{ftp_path}/{asm_name}_genomic.fna.gz", f"{asm_name}_genomic.fna.gz"),
                (f"{ftp_path}/{asm_name}_genomic.gff.gz", f"{asm_name}_genomic.gff.gz")
            ]

            for url, filename in files_to_download:
                filepath = os.path.join(outdir, filename)

                if file_exists(filepath):
                    print(f"{filename} already exists. Skipping.")
                    continue

                if not check_format(filename, os.path.splitext(filename)[1]):
                    logging.warning(f"Skipping {filename} due to format mismatch.")
                    continue
                    
                try:
                    with urllib.request.urlopen(url) as response:
                        if response.status != 200:
                            print(f"URL not accessible: {url}")
                            logging.warning(f"URL not accessible: {url}")
                            continue
                except Exception as url_error:
                    print(f"Error accessing URL {url}: {url_error}")
                    logging.warning(f"Error accessing URL {url}: {url_error}")
                    continue
                    
                try:
                    print(f"Downloading {filename} ...")
                    urllib.request.urlretrieve(url, filepath)
                    if file_exists(filepath):
                        print(f"Successfully downloaded {filename}")
                        logging.info(f"Downloaded {filename}")
                    else:
                        print(f"Failed to download {filename}")
                        logging.error(f"Failed to download {filename}")
                except Exception as download_error:
                    print(f"Error downloading {filename}: {download_error}")
                    logging.error(f"Error downloading {filename}: {download_error}")

        print(f"Reference genome download completed for {organism}")
        logging.info(f"Reference genome processing completed for {organism}")

    except Exception as e:
        logging.error(f"Error fetching reference genome: {e}")
        print(f"Error fetching reference genome: {e}")
# =================================================== #

# ==================== MAIN ==================== #
def my_parser() -> ap.Namespace:
    parser = ap.ArgumentParser(description="Pipeline to query and download BioProject data using Entrez.")
    parser.add_argument("-p", "--project", required=True, help="BioProject ID (e.g., PRJNA877658)")
    parser.add_argument("-e", "--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("-a", "--api_key", help="API key for NCBI Entrez (optional)")
    parser.add_argument("-o", "--outdir", default="bioproject_data", help="Output directory for all downloads")
    parser.add_argument("-i", "--info", action="store_true", help="Show BioProject general information")
    parser.add_argument("-d", "--databases", action="store_true", help="Show associated databases")
    parser.add_argument("-b", "--biosamples", action="store_true", help="Retrieve BioSamples information")
    parser.add_argument("-s", "--samples", action="store_true", help="Get SRA samples table")
    parser.add_argument("--download-all", action="store_true", help="Download all SRA FASTQ runs")
    parser.add_argument("--concat", action="store_true", help="Concatenate multiple FASTQs per sample")
    parser.add_argument("-r", "--reference", nargs='?', const="auto", metavar="ORGANISM",
                        help="Download reference genome. Without argument uses BioProject organism, or specify an organism")
    parser.add_argument("--max-ref-results", type=int, default=5,
                        help="Maximum number of reference results to show (default: 5)")
    return parser.parse_args()
# ============================================== #

# ==================== MAIN ==================== #
def main():
    args = my_parser()
    # Crear subdirectorios organizados
    sra_dir = os.path.join(args.outdir, "sra")
    concat_dir = os.path.join(args.outdir, "fastq_merged")
    biosample_dir = os.path.join(args.outdir, "biosamples")
    ref_dir = os.path.join(args.outdir, "reference_genome")
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(concat_dir, exist_ok=True)
    os.makedirs(biosample_dir, exist_ok=True)
    os.makedirs(ref_dir, exist_ok=True)

    try:
        configure_entrez(args.email, args.api_key)
    except ValueError as e:
        print(f"Configuration error: {e}")
        return

    print(f"\nStarting pipeline for BioProject: {args.project}")
    logging.info(f"Pipeline started for {args.project}")

    project_info = None
    if args.info:
        project_info = get_bioproject_info(args.project)
    if args.databases:
        get_associated_databases(args.project)
    if args.biosamples:
        get_biosample_info(args.project, output=os.path.join(biosample_dir, "biosamples.csv"))
    if args.samples:
        sra_df = get_sra_samples(args.project, output=os.path.join(sra_dir, "samples.csv"))
    else:
        sra_df = None

    if args.download_all or args.concat:
        if sra_df is None:
            sra_df = get_sra_samples(args.project, output=os.path.join(sra_dir, "samples.csv"))
        if sra_df is not None and (args.download_all or args.concat):
            download_and_concat_sra_grouped(sra_df, sra_dir, concat_dir)
        else:
            print("No SRA data available for download or concatenation")

    if args.reference:
        organism_to_search = args.reference
        if organism_to_search == "auto":
            if project_info is None:
                project_info = get_bioproject_info(args.project)
            if project_info:
                organism_to_search = project_info.get("Organism_Name")
                print(f"Using organism from BioProject: {organism_to_search}")
            else:
                print("Error: Could not retrieve BioProject info for organism")
                return
        if organism_to_search and organism_to_search != "N/A":
            get_reference_genome_info(organism_to_search, outdir=ref_dir,
                                      max_results=args.max_ref_results)
        else:
            print("Error: No organism specified or found in BioProject")

    logging.info("Pipeline execution completed successfully.")
    print(f"\nPipeline completed!")
    print(f"Log file saved as: {log_filename}")

if __name__ == "__main__":
    main()
