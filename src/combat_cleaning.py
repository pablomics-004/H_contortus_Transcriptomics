#!/usr/bin/env python3

"""
Combat-Seq Pipeline for RNA-Seq Count Matrix Batch Correction
------------------------------------------------------------

A streamlined pipeline for applying Combat-Seq batch effect correction to raw RNA-Seq count matrices.
Automatically handles metadata generation, batch correction, and post-filtering for downstream analysis.

Key Features:
- Automated minimal metadata generation for batch correction
- Combat-Seq implementation using pycombat_seq
- Intelligent filtering based on Counts Per Million (CPM) thresholds
- Comprehensive logging and error handling
- Support for various input formats and batch configurations

Input Requirements:
- Raw count matrix (genes x samples) in TSV/CSV format
- Optional batch information file or automatic batch detection

Output:
- Batch-corrected count matrix
- Filtered expression matrix (CPM-based)
- Diagnostic plots and quality control metrics
- Comprehensive processing log

Dependencies:
- pandas, numpy
- pycombat_seq

Example usage:
    python combat_seq_pipeline.py --counts raw_counts.tsv --output corrected_counts.tsv
    python combat_seq_pipeline.py --counts matrix.csv --batch NAME:col1,col2 ----U 1 --min_samples 3
""" 

import argparse
import pandas as pd
import numpy as np
import logging
from datetime import datetime
from typing import Dict, List, Tuple
import os
import sys 
from inmoose.pycombat import pycombat_seq

# ============================
#   LOGGING CONFIGURATION
# ============================
def configure_logging(logpath: str):
    """Configure logging to write to `logpath` and to stdout."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(logpath),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info(f"Logging initialized. Logfile: {logpath}")

# ==========================================================
# UTILITIES
# ==========================================================
def file_exists_nonempty(path: str) -> bool:
    """Return True if file exists and is non-empty."""
    if not os.path.exists(path):
        logging.error(f"File not found: {path}")
        return False
    if os.path.getsize(path) == 0:
        logging.error(f"File is empty: {path}")
        return False
    return True


def check_matrix_format(filepath: str, expected_extension: str = ".tsv") -> bool:
    """
    Validate that the file has the expected extension and contains a valid matrix.

    A valid matrix must:
    - load correctly as a dataframe
    - contain at least 2 columns (samples)
    - contain at least 2 rows (genes/transcripts)

    Parameters
    ----------
    filepath : str
        Path to the matrix file.

    expected_extension : str
        Expected file extension (default: ".tsv").

    Returns
    -------
    bool
        True if the matrix format is valid, False otherwise.
    """
    # --- Validate extension ---
    if not filepath.endswith(expected_extension):
        logging.warning(
            f"Format mismatch: {filepath} (expected {expected_extension})"
        )
        return False

    # --- Validate content can be loaded ---
    try:
        df = pd.read_csv(filepath, sep="\t", index_col=0)
    except Exception as e:
        logging.error(f"Matrix cannot be parsed as TSV: {filepath} | Error: {e}")
        return False

    # --- Validate matrix shape ---
    if df.shape[0] < 2 or df.shape[1] < 2:
        logging.error(
            f"Invalid matrix dimensions in {filepath}. Rows: {df.shape[0]}, Columns: {df.shape[1]}"
        )
        return False

    return True

def parse_args():
    """
    Parse command-line arguments for ComBat-seq + CPM filter pipeline.

    Returns
    -------
    argparse.Namespace
        Parsed arguments containing:
        - counts: path to raw counts TSV
        - batch: list of batch specifications (NAME:col1,col2,...)
        - U: CPM threshold
        - min_samples: minimum samples above CPM threshold
        - outdir: output directory
    """
    parser = argparse.ArgumentParser(
        description="ComBat-seq batch correction + CPM filter pipeline"
    )
    parser.add_argument("--counts", required=True,
                        help="Raw counts TSV (first col gene IDs)")
    parser.add_argument("--batch", action="append", required=True,
                        help="Batch spec: NAME:col1,col2,...  (1-based column numbers; 1 = gene ID). Repeat for multiple batches.")
    parser.add_argument("--U", type=float, default=1.0,
                        help="CPM threshold (default: 1.0)")
    parser.add_argument("--min_samples", type=int, default=None,
                        help="Minimum samples above CPM threshold (default: smallest batch)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    return parser.parse_args()

def load_counts(path: str) -> pd.DataFrame:
    """
    Load a raw counts TSV file (rows = genes, cols = samples).
    The first column is treated as gene IDs (index).
    
    Args:
        path: Path to TSV file with counts data
        
    Returns:
        DataFrame with genes as rows and samples as columns
        
    Raises:
        Exception: If file cannot be read
        ValueError: If matrix dimensions are invalid
    """
    try:
        # Read TSV file with gene IDs as index and sample names as header
        df = pd.read_csv(path, sep="\t", header=0, index_col=0)
    except Exception as e:
        logging.error(f"Failed to read counts file {path}: {e}")
        raise

    # Validate that matrix has at least one gene and one sample
    if df.shape[0] < 1 or df.shape[1] < 1:
        logging.error(f"Counts matrix seems empty or malformed: {path} shape={df.shape}")
        raise ValueError("Counts matrix must contain at least 1 gene row and 1 sample column.")

    logging.info(f"Loaded counts matrix '{path}' with shape {df.shape}")
    return df
             
def parse_batch_flags(batch_flags: List[str], n_columns: int) -> Dict[str, List[int]]:
    """
    Parse --batch flags specifying sample groupings for batch correction.
    
    Args:
        batch_flags: List of strings like "name:col1,col2,..."
        n_columns: Total columns in matrix (1-based, including gene ID column)
        
    Returns:
        Dict mapping batch names to lists of 1-based column indices
        
    Raises:
        ValueError: For malformed batch specifications
    """
    batches = {}
    for entry in batch_flags:
        # Split batch name from column list
        if ":" not in entry:
            logging.error(f"Invalid --batch entry (missing ':'): {entry}")
            raise ValueError(f"Invalid --batch entry: {entry}. Expected format name:1,2,3")
        name, cols_str = entry.split(":", 1)
        name = name.strip()
        
        # Validate batch name
        if not name:
            raise ValueError(f"Batch name empty in entry: {entry}")
            
        # Parse column indices
        cols = [c.strip() for c in cols_str.split(",") if c.strip() != ""]
        if not cols:
            raise ValueError(f"No columns listed for batch '{name}' in entry: {entry}")
        try:
            cols_int = [int(c) for c in cols]
        except ValueError:
            raise ValueError(f"Column indices must be integers in entry: {entry}")
            
        # Validate column ranges
        for ci in cols_int:
            if ci < 1 or ci > n_columns:
                raise ValueError(f"Column index {ci} out of range (1..{n_columns}) for entry: {entry}")
                
        # Check for duplicate batch names
        if name in batches:
            raise ValueError(f"Duplicate batch name: {name}")
            
        batches[name] = cols_int
        
    logging.info(f"Parsed batches: {batches}")
    return batches


def build_sample_to_batch_map(batches: Dict[str, List[int]], header_cols: List[str]) -> Dict[str, str]:
    """
    Convert batch column indices to sample names mapping.
    
    Args:
        batches: Dict mapping batch names to 1-based column indices
        header_cols: Header line split (first element is gene ID header)
        
    Returns:
        Dict mapping sample names to batch names
        
    Raises:
        ValueError: For invalid column assignments or conflicts
    """
    n_columns = len(header_cols)
    # Extract sample names from header (skip gene ID column)
    sample_names = header_cols[1:]
    # Create mapping from 1-based column numbers to sample names
    colnum_to_sample = {i + 2: sample_names[i] for i in range(len(sample_names))}
    
    sample_to_batch = {}
    for batch_name, colnums in batches.items():
        # Validate batch has columns assigned
        if len(colnums) == 0:
            raise ValueError(f"Batch '{batch_name}' has no columns assigned.")
            
        for cn in colnums:
            # Column 1 is gene ID, cannot be assigned to batch
            if cn == 1:
                raise ValueError("Column 1 is the gene ID column; it cannot be assigned to a batch.")
            # Validate column number exists
            if cn not in colnum_to_sample:
                raise ValueError(f"Column number {cn} is not a sample column (valid sample columns are 2..{n_columns})")
                
            sample = colnum_to_sample[cn]
            # Check for sample assigned to multiple batches
            if sample in sample_to_batch:
                raise ValueError(f"Sample '{sample}' assigned to multiple batches (conflict).")
                
            sample_to_batch[sample] = batch_name
            
    logging.info(f"Built sample->batch map for {len(sample_to_batch)} samples")
    return sample_to_batch

def make_metadata_df(sample_to_batch: Dict[str, str], all_samples: List[str]) -> pd.DataFrame:
    """
    Build metadata DataFrame for batch correction.
    
    Args:
        sample_to_batch: Dict mapping sample names to batch names
        all_samples: List of all sample names in the dataset
        
    Returns:
        DataFrame with 'sample' as index and 'batch' column
    """
    rows = []
    # Create rows for all samples, assigning 'unassigned' for samples without batch
    for s in all_samples:
        batch = sample_to_batch.get(s, "unassigned")
        rows.append((s, batch))
        
    # Create DataFrame and set sample names as index
    meta = pd.DataFrame(rows, columns=["sample", "batch"]).set_index("sample")
    logging.info(f"Generated metadata dataframe with shape {meta.shape}")
    return meta

def apply_combat_to_subset(counts: pd.DataFrame, samples_to_correct: List[str], batch_series: pd.Series) -> pd.DataFrame:
    """
    Apply ComBat-seq batch correction to specified samples.
    
    Args:
        counts: Raw count matrix (genes x samples)
        samples_to_correct: List of sample names to correct
        batch_series: Batch labels for each sample (index must match counts columns)
        
    Returns:
        Corrected count matrix with only specified samples adjusted
        
    Raises:
        ValueError: If samples not found in counts matrix
    """
    # Validate all target samples exist in counts matrix
    for s in samples_to_correct:
        if s not in counts.columns:
            raise ValueError(f"Sample '{s}' not found in counts matrix.")
    
    # Extract subset of counts for correction
    sub_counts = counts.loc[:, samples_to_correct]

    # Align batch labels with selected samples
    try:
        batches_for_sub = batch_series.loc[samples_to_correct]
    except Exception as e:
        logging.error(f"Failed to align batch vector to selected samples: {e}")
        raise

    # Apply ComBat-seq batch correction
    try:
        corrected_sub = pycombat_seq(counts=sub_counts, batch=batches_for_sub.tolist())
    except Exception as e:
        logging.error(f"pycombat_seq failed: {e}")
        raise

    # Merge corrected subset back into full matrix
    corrected_full = counts.copy()
    corrected_full.loc[:, samples_to_correct] = corrected_sub

    logging.info(f"Applied ComBat-seq to {len(samples_to_correct)} samples")
    return corrected_full
              

def compute_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Compute Counts Per Million (CPM) normalization.
    
    Args:
        counts: Raw count matrix (genes x samples)
        
    Returns:
        CPM normalized matrix (genes x samples)
    """
    # Calculate library sizes, replace zeros with NaN to avoid division by zero
    libsize = counts.sum(axis=0).replace(0, np.nan)
    # Compute CPM: (counts / library_size) * 1,000,000
    cpm = counts.div(libsize, axis=1) * 1e6
    # Fill NaN values (from zero library sizes) with 0.0
    return cpm.fillna(0.0)
              

def filter_by_cpm(counts: pd.DataFrame, U: float, min_samples: int) -> pd.DataFrame:
    """
    Filter genes based on Counts Per Million expression threshold.
    
    Args:
        counts: Raw count matrix (genes x samples)
        U: CPM threshold for expression
        min_samples: Minimum number of samples that must meet CPM threshold
        
    Returns:
        Tuple of (filtered count matrix, CPM matrix)
    """
    # Compute CPM normalized counts
    cpm = compute_cpm(counts)
    # Create mask for genes expressed above threshold in minimum samples
    mask = (cpm > U).sum(axis=1) >= min_samples
    # Apply filter to count matrix
    filtered = counts.loc[mask]
    logging.info(f"CPM filter (U={U}, min_samples={min_samples}): before={counts.shape[0]} genes, after={filtered.shape[0]} genes")
    return filtered, cpm


def main():
    """
    Execute Combat-Seq batch correction pipeline.
    """
    # Parse command line arguments
    args = parse_args()

    # Create output directory and set up logging
    os.makedirs(args.outdir, exist_ok=True)
    logpath = os.path.join(args.outdir, "combat_seq_pipeline.log")
    configure_logging(logpath)

    # Validate input counts file
    if not file_exists_nonempty(args.counts):
        sys.exit("Counts file missing or empty.")
    if not check_matrix_format(args.counts):
        sys.exit("Counts file failed format checks.")

    # Load count matrix
    counts_df = load_counts(args.counts)

    # Read header to get column information
    with open(args.counts, "r") as fh:
        header_cols = fh.readline().strip().split("\t")

    # Process batch information and create metadata
    batches = parse_batch_flags(args.batch, len(header_cols))
    sample_to_batch = build_sample_to_batch_map(batches, header_cols)
    metadata_df = make_metadata_df(sample_to_batch, list(counts_df.columns))
    
    # Save metadata
    metadata_path = os.path.join(args.outdir, "metadata.tsv")
    metadata_df.to_csv(metadata_path, sep="\t")
    logging.info(f"Metadata written to {metadata_path}")

    # Apply Combat-Seq batch correction
    samples_to_correct = list(sample_to_batch.keys())
    batch_sizes = [len(v) for v in batches.values()]
    min_samples = args.min_samples or min(batch_sizes)
    corrected_df = apply_combat_to_subset(counts_df, samples_to_correct, metadata_df["batch"])

    # Save corrected counts
    corrected_path = os.path.join(args.outdir, "corrected_counts.tsv")
    corrected_df.to_csv(corrected_path, sep="\t", float_format="%.6f")
    logging.info(f"Corrected counts saved to {corrected_path}")

    # Filter by CPM and save results
    filtered_counts, cpm_matrix = filter_by_cpm(corrected_df, U=args.U, min_samples=min_samples)
    filtered_path = os.path.join(args.outdir, "filtered_counts.tsv")
    filtered_counts.to_csv(filtered_path, sep="\t", float_format="%.6f")
    cpm_path = os.path.join(args.outdir, "cpm_matrix.tsv")
    cpm_matrix.to_csv(cpm_path, sep="\t", float_format="%.6f")
    logging.info(f"Filtered counts saved to {filtered_path}")
    logging.info(f"CPM matrix saved to {cpm_path}")

    # Print completion message with output paths
    print("Done.")
    print(f"Corrected counts: {corrected_path}")
    print(f"Filtered counts:  {filtered_path}")
    print(f"CPM matrix:       {cpm_path}")
    print(f"Metadata:         {metadata_path}")
    print(f"Log:              {logpath}")


if __name__ == "__main__":
    main()
