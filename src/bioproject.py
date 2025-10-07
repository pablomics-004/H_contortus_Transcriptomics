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
import re
import shlex
from urllib.request import urlretrieve, urlopen
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

def __to_srr(x: str) -> str:
    m = re.search(r'(SRR\d+)', x)
    return m.group(1) if m else str(x)

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
# =================================================== #

# ==================== BIOPROJECT ==================== #
def get_bioproject_info(bioproject_id: str) -> tuple[dict, int] | None:
    try:
        with Entrez.esearch(db="bioproject", term=bioproject_id) as handle:
            record = Entrez.read(handle)

        if not record["IdList"]:
            raise ValueError(f"Couldn't find BioProject {bioproject_id}")

        uid = record["IdList"][0]

        with Entrez.esummary(db="bioproject", id=uid) as summary:
            results = Entrez.read(summary)
        result = results["DocumentSummarySet"]["DocumentSummary"][0]

        msg = (
            f"\n=== BioProject Information ===\n"
            f"Title: {result.get('Project_Title', 'N/A')}\n"
            f"Accession: {result.get('Project_Acc', 'N/A')}\n"
            f"Description: {result.get('Project_Description', 'N/A')}\n"
            f"Organism: {result.get('Organism_Name', 'N/A')}\n"
            f"Data Type: {result.get('Project_Data_Type', 'N/A')}\n"
            f"Submitter: {result.get('Submitter_Organization', 'N/A')}\n"
            f"Registration Date: {result.get('Registration_Date', 'N/A')}\n"
            f"URL: https://www.ncbi.nlm.nih.gov/bioproject/{bioproject_id}\n"
        )
        print(msg)
        logging.info(msg)
        return result, uid

    except Exception as e:
        logging.error(f"Error fetching BioProject info: {e}")
        print(f"Error fetching BioProject info: {e}")
        return None
# =================================================== #

# ==================== DATABASE LINKS ==================== #
def get_associated_databases(uid: str) -> list[str] | None:
    try:
        with Entrez.elink(dbfrom="bioproject", id = uid) as handle:
            records = Entrez.read(handle)

        linked_dbs = sorted({
            db['DbTo']
            for linkset in records if 'LinkSetDb' in linkset
            for db in linkset['LinkSetDb']
        })

        if linked_dbs:
            msg = "=== Associated Databases ===\n" + "\n".join(f"- {db}" for db in linked_dbs)
        else:
                msg = f'No associated data fot {uid}'
        
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
            record = Entrez.read(handle)
        return {
            'BioSample_ID': bid,
            'Accession': summary.get('Accession', 'N/A'),
            'Organism': summary.get('Organism', 'N/A'),
            'Title': summary.get('Title', 'N/A'),
            'Attributes': summary.get('Attributes', 'N/A'),
            'SubmissionDate': summary.get('SubmissionDate', 'N/A')
        }
    except Exception as inner_e:
        logging.warning(f"Error retrieving BioSample {bid}: {inner_e}")
        return
# =================================================== #

# ==================== BIOSAMPLES ==================== #
def get_biosample_info(uid: str, output: str) -> pd.DataFrame | None:
    try:
        with Entrez.elink(dbfrom='bioproject', id=uid, db='biosample') as handle:
            records = Entrez.read(handle)
        
        biosample_ids = [
            link['Id']
            for linkset in records
            for db in linkset.get('LinkSetDb', [])
            for link in db.get('Link', [])
            if 'Id' in link
        ]

        if not biosample_ids:
            print("No BioSamples linked to this BioProject.")
            return
        
        print(f"{len(biosample_ids)} BioSamples found. Downloading metadata...")

        results = [data for bid in biosample_ids if (data := fetch_biosamples(bid))]

        if results:
            df = pd.DataFrame(results)
            df.to_csv(output, index=False)
            print(f"Retrieved {len(df)} BioSamples. Saved to {output}")
            logging.info(f"Retrieved {len(df)} BioSamples.")
            return
        
        print('No valid BioSamples were retrieved.')
    
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
        
        # Extracting the first run Accession within the runs list
        runs_text = summary.get('Runs', dft:='N/A')
        run_accession = runs_text.split(',')[0] if runs_text != dft else dft
        return {
            'RunAccession' : run_accession,
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
# =============================================== #

# ==================== SRA ==================== #
def get_sra_samples(bioproject_id: str, output: str) -> pd.DataFrame | None:
    try:
        with Entrez.elink(dbfrom="bioproject", db="sra", id=bioproject_id) as handle:
            records = Entrez.read(handle)
        sra_ids = {
            link['Id']
            for linkset in records
            for db in linkset.get('LinkSetDb', [])
            for link in db.get('Link', [])
            if 'Id' in link
        }

        if not sra_ids:
            msg = f"No SRA entries linked to BioProject {bioproject_id}"
            print(msg); logging.info(msg)
            return
        
        msg = f"Found {len(sra_ids)} SRA entries linked to {bioproject_id}"
        print(msg); logging.info(msg)

        sample_data = [data for sra_id in sra_ids if (data := fetch_sra(sra_id))]

        df = pd.DataFrame(sample_data)
        df.to_csv(output, index=False)
        msg = f"Saved SRA sample info to {output} ({len(df)} valid runs)"
        print(msg); logging.info(msg)
        return

    except Exception as e:
        logging.error(f"Error fetching SRA samples for {bioproject_id}: {e}")
        print(f"Error fetching SRA samples: {e}")
        return
# =================================================== #

# ==================== SRA DOWNLOAD & CONCAT ==================== #
def download_and_concat_sra_grouped(sra_df: pd.DataFrame, sra_dir: str, concat_dir: str, bash: bool = False) -> None:

    os.makedirs(sra_dir, exist_ok=True)
    os.makedirs(concat_dir, exist_ok=True)

    grouped = defaultdict(list)
    layout_map = {}

    for _, row in sra_df.iterrows():
        grouped[row['Title']].append(row['RunAccession'])
        layout_map[row['Title']] = row['LibraryLayout']

    if bash: 
        here = os.path.dirname(os.path.abspath(__file__))
        script1 = os.path.join(here, 'pf_fastq.sh')
        script2 = os.path.join(here, 'concatenate_fastq.sh')

        # Validation for my own
        for s in (script1, script2):
            if not os.path.isfile(s):
                raise FileNotFoundError(f"Required bash script not found: {s}")
            if not os.access(s, os.X_OK):
                raise PermissionError(f"Bash script is not executable: {s}")

        manifest_path = os.path.join(sra_dir, f'manifest_runs.{int(time.time())}.tsv')

        for sample, runs in grouped.items():
            run_ids = [__to_srr(r) for r in runs]
            layout = str(layout_map.get(sample, 'SINGLE')).upper()
            paired = (layout == 'PAIRED')

            cmmd = [
                script1,
                '-s', sra_dir,
                '-c', concat_dir,
                '-S', sample,
                '-o', sample,
                '-m', manifest_path,
            ]
            if paired:
                cmmd.append('-p')
            cmmd += run_ids

            print("[BASH] ", " ".join(shlex.quote(x) for x in cmmd)); subprocess.run(cmmd, check=True)

        # Add -c if you want to use fast_cat
        cmmd2 = [script2, '-m', manifest_path,'-f']
        print("[BASH] ", " ".join(shlex.quote(x) for x in cmmd2)); subprocess.run(cmmd2, check=True)

        print(f"[INFO] Concatenation done from manifest: {manifest_path}")
        return

    for sample_name, runs in grouped.items():
        r1_files, r2_files = [], []
        
        for run in runs:
            run_id = __to_srr(run)

            run_dir = os.path.join(sra_dir, run_id)
            os.makedirs(run_dir, exist_ok=True)

            r1_path = os.path.join(run_dir, f"{run_id}_1.fastq.gz")
            r2_path = os.path.join(run_dir, f"{run_id}_2.fastq.gz") if layout_map[sample_name] == "PAIRED" else None

            if not file_exists(r1_path) or (r2_path and not file_exists(r2_path)):
                print(f"Prefetching {run_id} ...")
                subprocess.run(["prefetch", "--output-directory", sra_dir, run_id], check=True)
                print(f"Converting {run_id} to FASTQ.gz ...")
                subprocess.run(["fastq-dump", run_id, "-O", run_dir, "--split-files", "--gzip"], check=True)
            else:
                print(f"{run_id} already downloaded. Skipping.")

            if file_exists(r1_path) and check_format(r1_path, ".fastq.gz"):
                r1_files.append(r1_path)
            if r2_path and file_exists(r2_path) and check_format(r2_path, ".fastq.gz"):
                r2_files.append(r2_path)

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

            urls = [f"{ftp_path}/{asm_name}_genomic.fna.gz",f"{ftp_path}/{asm_name}_genomic.gff.gz"]
            filenames = [f"{asm_name}_genomic.fna.gz",f"{asm_name}_genomic.gff.gz"]

            for url, filename in zip(urls, filenames):
                filepath = os.path.join(outdir, filename)

                if file_exists(filepath):
                    print(f"{filename} already exists. Skipping.")
                    continue

                if not check_format(filename, os.path.splitext(filename)[1]):
                    logging.warning(f"Skipping {filename} due to format mismatch.")
                    continue
                    
                try:
                    with urlopen(url) as response:
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
                    urlretrieve(url, filepath)
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
    parser.add_argument("-p", "--project", type=str, required=True, help="BioProject ID (e.g., PRJNA877658)")
    parser.add_argument("-e", "--email", type=str, required=True, help="Email for NCBI Entrez")
    parser.add_argument("-a", "--api_key", type=str, help="API key for NCBI Entrez (optional)")
    parser.add_argument("-o", "--outdir", type=str, default="bioproject_data", help="Output directory for all downloads")
    parser.add_argument("-i", "--info", action="store_true", help="Show BioProject general information")
    parser.add_argument("-d", "--databases", action="store_true", help="Show associated databases")
    parser.add_argument("-b", "--biosamples", action="store_true", help="Retrieve BioSamples information")
    parser.add_argument("-s", "--samples", action="store_true", help="Get SRA samples table")
    parser.add_argument("--download_all", action="store_true", help="Download all SRA FASTQ runs")
    parser.add_argument("--concat", action="store_true", help="Concatenate multiple FASTQs per sample")
    parser.add_argument("-u", "--use_bash", action="store_true", help="Allows the downloading and concatenation of sra by a bash script")
    parser.add_argument("-r", "--reference", nargs='?', const="auto", metavar="ORGANISM",
                        help="Download reference genome. Without argument uses BioProject organism, or specify an organism")
    parser.add_argument("--max_ref_results", type=int, default=5,
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
            download_and_concat_sra_grouped(sra_df, sra_dir, concat_dir, args.use_bash)
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
