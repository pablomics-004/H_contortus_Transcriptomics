# Comparative transcriptomics of Haemonchus contortus: ivermectin susceptibility and resistance (IVM-S/IVM-R)
Date: 29/08/2025

This repository contains the bioinformatics pipeline and transcriptomic analysis of Haemonchus contortus, focusing on the molecular basis of ivermectin resistance.
The project compares ivermectin-susceptible (IVM-S) and ivermectin-resistant (IVM-R) strains through RNA-seq differential expression analysis. This project reproduces and extends the analysis from [Reyes-Guerrero _et.al_ (2023)](https://doi.org/10.3390/pathogens12030499).

**Participants**:  
- Pablo Salazar-Mendez <pablosm@lcg.unam.mz>
- Ashley Yael Montiel-Vargas <yaelmont@lcg.unam.mx>

___
## 📌 Project Description
This project performs a **comparative transcriptomic analysis of Haemonchus contortus**, a hematophagous gastrointestinal nematode that causes severe economic losses in livestock.
The study focuses on comparing **ivermectin-susceptible (IVM-S)** and **ivermectin-resistant (IVM-R)** strains to identify **differentially expressed genes (DEGs)** that may underlie anthelmintic resistance mechanisms. Understanding these molecular signatures is crucial for improving parasite control strategies and delaying the spread of drug resistance in field populations.

Ultimately, this project aims to provide a **reproducible computational framework for transcriptomic studies** , enabling future research on drug resistance, gene regulation, and parasite adaptation mechanisms.
___

## 🎯 Objectives

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
3. Quality control and trimming of raw reads (FastQC, Trimmomatic, Cutadpat).  
4. Mapping and quantification of transcripts (HISAT2, Salmon/Kallisto, BWA).  
5. Differential expression analysis (DESeq2, EdgeR).  
6. Functional annotation (BLAST, InterProScan, GO/KEGG).  
7. Integration and interpretation of results.  

___ 

## Expected Results
 
- Lists of differentially expressed genes between susceptible and resistant isolates.  
- Functional insights into molecular pathways associated with ivermectin resistance.  
- Reproducible pipelines and scripts deposited in GitHub.  
- A documented workflow for future studies on helminth transcriptomics.  

___

## Requirements Specification

**Functional Requirements**  

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
Step 3: Map_and_Quantify(assembly, reads)
Step 4: Differential_Expression(quant_data)
Step 5: Functional_Annotation(diff_exp_genes)
Step 6: Generate_Reports(results)
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
- **Description**: The researcher provides raw RNA-seq data for IVM-S and IVM-R strains. The pipeline validates and processes the data, performs quantification, differential expression, and annotation, and outputs interpretable results.  
- **Main Flow**:  
  1. Input RNA-seq data.  
  2. Quality check and preprocessing.    
  3. Differential expression analysis.  
  4. Functional annotation.  
  5. Results reporting.  

- **Alternative Flows**:  
  - If input data is corrupted → System halts and reports error.  
  - If assembly fails → Provide error log and suggest parameter adjustments.  
  - If no differential expression is detected → Report and flag possible biological/technical causes.  

___

 ## ⚙️ Installation and Usage

Instructions for environment setup and execution will be provided as the project develops.
Due to the large file sizes of RNA-seq datasets (>2 GB), raw data are not stored directly in this repository but are accessible via NCBI.
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
## 🧾 License

This project is distributed under the MIT License.
See the LICENSE file for details.
___

## 🧩 Citation

If you use this work in your research, please cite:

Salazar-Mendez P, Montiel-Vargas AY. (2025). Comparative transcriptomics of Haemonchus contortus: ivermectin susceptibility and resistance (IVM-S/IVM-R) . [GitHub Repository]. https://github.com/pablomics-004/haemonchus_transcriptome_analysis
___

## 📬 Contact

For questions, suggestions, or issues, please open an Issue in this repository or contact:

Pablo Salazar-Mendez — <pablosm@lcg.unam.mx>

Ashley Yael Montiel-Vargas — <yaelmont@lcg.unam.mx>

___ 

## 📚 References

Reyes-Guerrero, D. E., Jiménez-Jacinto, V., Alonso-Morales, R. A., Alonso-Díaz, M. Á., Maza-Lopez, J., Camas-Pereyra, R., Olmedo-Juárez, A., Higuera-Piedrahita, R. I., & López-Arellano, M. E. (2023). Assembly and Analysis of Haemonchus contortus Transcriptome as a Tool for the Knowledge of Ivermectin Resistance Mechanisms. Pathogens, 12(3), 499. https://doi.org/10.3390/pathogens12030499 

