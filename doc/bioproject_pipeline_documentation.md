# 🧬 BioProject Entrez Data Retrieval and SRA Processing Pipeline

## Overview

This script automates querying and downloading data associated with an NCBI BioProject. It retrieves metadata (BioProject, BioSamples, and SRA runs), downloads sequencing data (FASTQ files), and optionally fetches the reference genome.

It integrates multiple steps into a single, reproducible workflow:

- Query BioProject metadata from NCBI Entrez.
- Retrieve associated BioSamples and SRA runs.
- Download FASTQ files using **prefetch** and **fastq-dump**.
- Concatenate per-sample reads into unified FASTQ files.
- Optionally, fetch the reference genome.
- All actions are logged with timestamps in a detailed log file.

---

## Requirements

### Dependencies

- **Python ≥ 3.9**
- **Biopython**
- **pandas**
- **NCBI SRA Toolkit** (`prefetch`, `fastq-dump`)
- `gzip` / `zcat` utilities (Linux/macOS)
- **sra-toolkit**

## FILE ORGANIZATION

``` bash
bioproject_data/
├── biosamples/               # BioSample metadata (CSV)
├── sra/                      # Individual downloaded SRA runs
├── fastq_merged/             # Concatenated FASTQ files per sample
├── reference_genome/         # Reference genome sequences and annotations
└── bioproject_pipeline_YYYYMMDD_HHMMSS.log
```

## KEY COMPONENTS

1. Logging Configuration

    - Creates console and file logs for reproducibility and debugging.

2. Utility Functions

    - ile_exists(): Checks if a file exists and is non-empty.

    - check_format(): Verifies expected file extensions.

    - concat_gz(): Concatenates .fastq.gz files efficiently using zcat and gzip.

3. Entrez Configuration
    - Sets up NCBI Entrez with a valid email (required) and optional API key.
```python
configure_entrez(email: str, api_key: Optional[str])
```

4. BioProject Metadata

```python
get_bioproject_info(bioproject_id: str)
```
Retrieves key BioProject information:

- Title

- Accession

- Description

- Data type

- Submitter organization

- Registration date

- Returns both the metadata and UID for subsequent queries.

5. Linked Databases
```python
get_associated_databases(uid: str)
```

Lists all databases associated with a BioProject (e.g., BioSample, SRA, Assembly).

6. BioSample Metadata

```python
fetch_biosamples(bioproject_uid: str, output: str = None)
```

Retrieves BioSample metadata:

- Accession ID

- Organism name

- Title and description

- Submission date

- SRA identifiers

- Attributes

Saves results as a CSV and returns a pandas DataFrame.

7. SRA Retrieval

```python
get_sra_by_bioproject(bioproject_uid: str, output_csv: str)
```

- Fetches all SRA run accessions (SRR) linked to each BioSample, including sequencing layout (SINGLE/PAIRED).

Saved as sra_by_biosample.csv with columns:

- BioSample

- SRR_List (comma-separated)

- ibraryLayout

8. Reference Genome

```python
get_reference_genome_info(organism: str, outdir="reference_genome", max_results=5)
```

Searches NCBI Assembly for the latest reference genomes of the target organism.

Downloads:

- *_genomic.fna.gz (FASTA sequence)

- *_genomic.gff.gz (annotations)

9. SRA Download and Concatenation

```python
download_and_concat_sra_grouped(sra_df: DataFrame, sra_dir: str, concat_dir: str)
```

Downloads SRR runs using the SRA Toolkit (prefetch and fastq-dump), decompresses into .fastq.gz, and merges per sample:

- Paired-end: _1.fastq.gz and _2.fastq.gz concatenated separately.

- Single-end: concatenates single read files.

All actions are logged.

10. Main Function (CLI)

Handles all command-line options via argparse and orchestrates the workflow.


## Command-Line Usage

### Basic Structure

```bash
python bioproject_pipeline.py -p PRJXXXXXXX -e your_email@example.com [OPTIONS]
```


#### Retrieve BioProject info only

```bash
python bioproject_pipeline.py -p PRJNA877658 -e your_email@example.com -i
```

#### Retrieve BioSamples and SRA runs:
```bash
python bioproject_pipeline.py -p PRJNA877658 -e your_email@example.com -b -s

```

#### Download all SRA runs and concatenate per sample:

```bash
python bioproject_pipeline.py -p PRJNA877658 -e your_email@example.com --download-all --concat

```

#### Download BioSamples, SRA runs, and the reference genome automatically

```bash
python bioproject_pipeline.py -p PRJNA877658 -e your_email@example.com -b -s -r auto

```

## Outputs

- **biosamples.csv** — Table of BioSample metadata.

- **sra_by_biosample.csv** — Table linking BioSamples to SRR runs and library layout.

- **FASTQ files** — One per SRR run in sra/, merged by sample in fastq_merged/.

- **Reference genome** files — Downloaded FASTA and GFF files in reference_genome/.

- **Log file** — Records all actions, warnings, and errors with timestamps.

## NOTES

- Always use a valid email for NCBI Entrez.

- For high-volume queries, register and use an NCBI API key.

- Ensure the SRA Toolkit commands (prefetch, fastq-dump) are in your system PATH.

- If encountering connection errors, add short delays between Entrez queries to respect NCBI rate limits.