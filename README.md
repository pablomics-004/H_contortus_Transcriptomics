# Transcriptome Assembly and Analysis of Haemonchus contortus (IVM-S/IVM-R)

Date: 29/08/2025

**Participants**:  
- Pablo Salazar-Mendez <pablosm@lcg.unam.mz>
- Ashley Yael Montiel-Vargas <yaelmont@lcg.unam.mx>

___
## 📌 Project Description
This project aims to perform the *transcriptome assembly and analysis* of Haemonchus contortus, a gastrointestinal parasite of veterinary importance.  
The goal is to compare *ivermectin-susceptible (IVM-S)* and *ivermectin-resistant (IVM-R)* strains to identify differences in gene expression associated with drug resistance.
___

## 🎯 Objectives
- Assemble and assess the quality of the H. contortus transcriptome.  
- Compare gene expression profiles between IVM-S and IVM-R strains.  
- Identify differentially expressed genes related to drug resistance.  
- Generate reproducible resources for downstream analyses.

___

## 🧬 Data
- Organism: Haemonchus contortus  
- Conditions:  
  - Ivermectin-susceptible strain (IVM-S)
    Transcriptome availabilty (SRR files):
    - [SRR21518936](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518936&display=download) (`1.9Gb`)
    - [SRR21518937](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518937&display=download) (`1.6Gb`)
    - [SRR21518938](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518938&display=download) (`1.6Gb`)
  - Ivermectin-resistant strain (IVM-R)
    Transcriptome availabilty (SRR files):
    - [SRR21518939](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518939&display=download) (`2Gb`)
    - [SRR21518940](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518940&display=download) (`2.7Gb`)
    - [SRR21518941](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR21518941&display=download) (`1.5Gb`)
- Data type: RNA-seq (Illumina)
___

## Work Plan / Timeline

You can check all the workplan and progress on the project associated with this repository <https://github.com/users/pablomics-004/projects/2>
___

## Methodology

**Research Questions**  
- What transcriptional differences exist between IVM-S and IVM-R *H. contortus* isolates?  
- Which genes or pathways are associated with ivermectin resistance?  
- How can transcriptome data support better understanding of drug resistance mechanisms in parasitic nematodes?  

**Steps**  
1. Data source identification (RNA-seq Illumina data from IVM-S and IVM-R strains).  
2. Download and preparation of raw sequencing data.  
3. Quality control and trimming of raw reads (FastQC, Trimmomatic).  
4. De novo transcriptome assembly (Trinity).  
5. Mapping and quantification of transcripts (HISAT2, Salmon/Kallisto).  
6. Differential expression analysis (DESeq2, EdgeR).  
7. Functional annotation (BLAST, InterProScan, GO/KEGG).  
8. Integration and interpretation of results.  

___ 

## Expected Results

- A high-quality transcriptome assembly of *H. contortus* (IVM-S and IVM-R).  
- Lists of differentially expressed genes between susceptible and resistant isolates.  
- Functional insights into molecular pathways associated with ivermectin resistance.  
- Reproducible pipelines and scripts deposited in GitHub.  
- A documented workflow for future studies on helminth transcriptomics.  

___

## Requirements Specification

**Functional Requirements**  
- The pipeline must assemble transcriptome data from raw RNA-seq reads.  
- It must quantify transcript expression and perform differential expression analysis.  
- It must annotate genes and provide biological context for differentially expressed transcripts.  
- It must generate reproducible scripts for all steps of the workflow.  

**Non-Functional Requirements**  
- Scripts must be written in reproducible languages (R, Python, Bash).  
- The workflow should be modular and easy to adapt to other transcriptomic datasets.  
- Output should be well-documented and structured for downstream analyses.  
- Performance should handle large RNA-seq datasets efficiently.  

___

## Analysis and Design

The project will integrate a modular pipeline combining existing bioinformatics tools. Workflow management will focus on data validation, error handling, and reproducibility.  

**Pipeline Design (simplified pseudocode):**
Main Workflow (Transcriptome_Analysis):
Input: Raw RNA-seq reads (IVM-S, IVM-R)
Step 1: Quality_Control(reads)
Step 2: Trim_Reads(reads)
Step 3: Assemble_Transcriptome(trimmed_reads)
Step 4: Map_and_Quantify(assembly, reads)
Step 5: Differential_Expression(quant_data)
Step 6: Functional_Annotation(diff_exp_genes)
Step 7: Generate_Reports(results)
Output: Expression profiles, annotation, reproducible scripts

**Use Case: Transcriptome Analysis**
     +---------------+
     |   Researcher  |
     +-------+-------+
             |
             | 1. Provides RNA-seq data
             v
     +-------+-------+
     | Transcriptome |
     |  Pipeline     |
     +---------------+

- **Actor**: Researcher  
- **Description**: The researcher provides raw RNA-seq data for IVM-S and IVM-R strains. The pipeline validates and processes the data, performs assembly, quantification, differential expression, and annotation, and outputs interpretable results.  
- **Main Flow**:  
  1. Input RNA-seq data.  
  2. Quality check and preprocessing.  
  3. Transcriptome assembly and quantification.  
  4. Differential expression analysis.  
  5. Functional annotation.  
  6. Results reporting.  

- **Alternative Flows**:  
  - If input data is corrupted → System halts and reports error.  
  - If assembly fails → Provide error log and suggest parameter adjustments.  
  - If no differential expression is detected → Report and flag possible biological/technical causes.  

___

## 📂 Repository Structure

```
├── data/       # Raw or processed data 
├── doc/        # Documentation (project notes, references, reports)
├── results/    # Processed results (tables, figures, reports)
├── src/        # Source code and analysis scripts (bash, Python, R)
├── lib/        # External libraries, modules, or custom functions
├── test/       # Test scripts and validation files
├── tmp/        # Temporary or intermediate files
├── LICENSE     # License information for the project
└── README.md   # Project description and instructions
```
___