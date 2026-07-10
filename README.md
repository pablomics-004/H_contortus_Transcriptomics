# Comparative Transcriptomics of _Haemonchus contortus_: Ivermectin Susceptibility and Resistance (IVM-S/IVM-R)

[![Python](https://img.shields.io/badge/Python-3.x-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Bash](https://img.shields.io/badge/Bash-Scripting-4EAA25?logo=gnu-bash&logoColor=white)]()
[![PyDESeq2](https://img.shields.io/badge/PyDESeq2-Differential_Expression-8A2BE2)](https://pydeseq2.readthedocs.io/)
[![Biopython](https://img.shields.io/badge/Biopython-(Bio.Entrez)-F3C623?logo=python&logoColor=3776AB)](https://biopython.org/)
[![Data: RNA-seq](https://img.shields.io/badge/Data-RNA--seq-lightgrey.svg)]()

**Date:** August 29th, 2025

This repository contains the bioinformatics pipeline and transcriptomic analysis of *Haemonchus contortus*, focusing on the molecular basis of ivermectin resistance. The project compares ivermectin-susceptible (IVM-S) and ivermectin-resistant (IVM-R) strains through RNA-seq differential expression analysis. This project reproduces and extends the analysis from [1].

**Participants:** 

* Pablo Salazar-Mendez <pablosm@lcg.unam.mx>
* Ashley Yael Montiel-Vargas <yaelmont@lcg.unam.mx>

---

## Project Description

This project performs a comparative transcriptomic analysis of *Haemonchus contortus*, a hematophagous gastrointestinal nematode that causes severe economic losses in the livestock industry [1]. 

The study focuses on comparing ivermectin-susceptible (IVM-S) and ivermectin-resistant (IVM-R) strains to identify differentially expressed genes (DEGs) that may underlie anthelmintic resistance mechanisms. Understanding these molecular signatures is crucial for improving parasite control strategies and delaying the spread of drug resistance in field populations [1].

Ultimately, this project aims to provide a reproducible computational framework for transcriptomic studies, utilizing robust Python libraries (such as PyDESeq2 for statistical modeling and Bio.Entrez for programmatic data retrieval), enabling future research on drug resistance, gene regulation, and parasite adaptation mechanisms [1].

---

## Objectives

* Compare global gene expression profiles between IVM-S and IVM-R strains.  
* Identify specific differentially expressed genes related to drug resistance.  
* Generate reproducible computational resources and pipelines for downstream analyses.

---

## Data Specifications

* **Organism:** *Haemonchus contortus* * **Data Type:** RNA-seq (Illumina)
* **Conditions & Availability:** * **Ivermectin-susceptible strain (IVM-S):**
    * [SRR21518936](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518936&display=download) (`1.9 Gb`)
    * [SRR21518937](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518937&display=download) (`1.6 Gb`)
    * [SRR21518938](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518938&display=download) (`1.6 Gb`)
  * **Ivermectin-resistant strain (IVM-R):**
    * [SRR21518939](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518939&display=download) (`2.0 Gb`)
    * [SRR21518940](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518940&display=download) (`2.7 Gb`)
    * [SRR21518941](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518941&display=download) (`1.5 Gb`)

---

## Work Plan and Timeline

All project tracking, milestones, and task progression can be monitored via the associated GitHub Project board:  
[Project Management Dashboard](https://github.com/users/pablomics-004/projects/2)

---

## Methodology

### Research Questions  
* What transcriptional differences exist between IVM-S and IVM-R *H. contortus* isolates?  
* Which genetic pathways are directly associated with ivermectin resistance?  
* How can whole-transcriptome data support a better understanding of drug resistance mechanisms in parasitic nematodes?  

### Analytical Steps  
1. **Data Acquisition:** Identification and programmatic download of raw RNA-seq Illumina data via `Bio.Entrez` and SRA-Toolkit.  
2. **Preprocessing:** Quality control and adapter trimming of raw reads.  
3. **Alignment & Quantification:** Mapping reads to the reference genome and quantifying transcripts.  
4. **Statistical Analysis:** Differential expression analysis utilizing Python-native models (`PyDESeq2`).  
5. **Functional Annotation:** Assigning biological context via BLAST, InterProScan, and GO/KEGG pathway enrichment.  
6. **Synthesis:** Integration and biological interpretation of the computational results.  

--- 

## Expected Results
 
* Comprehensive lists of differentially expressed genes between susceptible and resistant isolates.  
* Functional insights into the molecular pathways driving ivermectin resistance.  
* Highly reproducible analytical pipelines and scripts hosted in this repository.  
* A strictly documented workflow establishing a baseline for future helminth transcriptomic studies.  

---

## Requirements Specification

### Functional Requirements  
* The system must accurately quantify transcript expression and perform rigorous differential expression analysis natively in Python (`PyDESeq2`).  
* The pipeline must functionally annotate genes to provide biological context for the identified transcripts.  
* The workflow must automatically generate reproducible scripts and logs for all processing steps.  

### Non-Functional Requirements  
* Codebases must be written in standard, reproducible languages (Python, R, Bash).  
* The workflow architecture must be modular, allowing for easy adaptation to subsequent transcriptomic datasets.  
* Output files and directories must be well-documented and logically structured for seamless downstream integration.  
* Algorithmic performance must be optimized to handle large-scale RNA-seq datasets efficiently.  

---

## Analysis and Design

The project integrates a modular pipeline combining standard bioinformatics tools. The workflow management heavily emphasizes data validation, error handling, and computational reproducibility.  

**Pipeline Architecture (Simplified Pseudocode):**
```text
Main Workflow (Transcriptome_Analysis):
Input: Raw RNA-seq reads (IVM-S, IVM-R)
Step 1: Quality_Control(reads)
Step 2: Trim_Reads(reads)
Step 3: Map_and_Quantify(assembly, reads)
Step 4: Differential_Expression(quant_data) # Handled via PyDESeq2
Step 5: Functional_Annotation(diff_exp_genes)
Step 6: Generate_Reports(results)
Output: Expression profiles, functional annotations, reproducible scripts
```

## System Use Case: Transcriptome Pipeline

```text
     +---------------+
     |  Researcher   |
     +-------+-------+
             |
             | 1. Provides raw RNA-seq data
             v
     +-------+-------+
     | Transcriptome |
     |   Pipeline    |
     +---------------+
```

* Actor: Researcher

* Description: The researcher inputs raw RNA-seq data for both strains. The automated pipeline validates the data, processes it through quantification and differential expression, applies functional annotation, and outputs interpretable biological results.

* Alternative Flows:
  * Data Corruption: System halts execution and generates an error report.
  * Alignment Failure: System outputs an error log and prompts for parameter adjustments.
  * Null Differential Expression: System reports the lack of significant DEGs and flags potential technical or biological variances.
 
## Installation and Usage

Detailed instructions for environment configuration and pipeline execution will be updated as the codebase matures.

>Note: Due to the substantial file sizes characteristic of raw RNA-seq datasets (`>2 GB` per file), .fastq files are not hosted directly within this repository. They must be downloaded via the NCBI links provided in the Data section.

## Repository Structure

```text
.
├── data
├── doc
│   ├── Abstract_RNAgels_2.md
│   ├── bioproject_pipeline_documentation.md
│   └── RNAngels_final_report.md
├── LICENSE
├── README.md
├── results
│   ├── annotation_enrichment
│   ├── combat
│   ├── filtered_annotation
│   └── Images
├── src
│   ├── bioproject_pipeline.py
│   ├── combat_cleaning.py
│   ├── concat_gz.sh
│   ├── differential_expression_analysis.py
│   ├── enrichment.py
│   ├── fastqc_precheck.sh
│   ├── process_trionate.sh
│   ├── __pycache__
│   │   └── visualization.cpython-312.pyc
│   ├── rna_cuantification.sh
│   ├── run_prefetch.sh
│   ├── temp_script.sh
│   ├── testing_mapper.ipynb
│   ├── trimming.sh
│   └── visualization.py
└── tmp
```

## Citation

If you utilize this workflow or these findings in your research, please cite the repository using the following BibTeX format:

```text
@misc{salazar_montiel_2025_haemonchus,
  author       = {Salazar-Mendez, Pablo and Montiel-Vargas, Ashley Yael},
  title        = {Comparative transcriptomics of Haemonchus contortus: ivermectin susceptibility and resistance (IVM-S/IVM-R)},
  year         = {2025},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{[https://github.com/pablomics-004/haemonchus_transcriptome_analysis](https://github.com/pablomics-004/haemonchus_transcriptome_analysis)}}
}
```

## Contact Information

For inquiries, methodological suggestions, or technical issues, please open an Issue within this repository or contact the authors directly:

* Pablo Salazar-Mendez: <pablosm@lcg.unam.mx>
* Ashley Yael Montiel-Vargas: <yaelmont@lcg.unam.mx>

## References

[1]  Reyes-Guerrero, D. E., Jiménez-Jacinto, V., Alonso-Morales, R. A., Alonso-Díaz, M. Á., Maza-Lopez, J., Camas-Pereyra, R., Olmedo-Juárez, A., Higuera-Piedrahita, R. I., & López-Arellano, M. E. (2023). Assembly and Analysis of Haemonchus contortus Transcriptome as a Tool for the Knowledge of Ivermectin Resistance Mechanisms. Pathogens, 12(3), 499. [https://doi.org/10.3390/pathogens12030499](https://www.mdpi.com/2076-0817/12/3/499)
