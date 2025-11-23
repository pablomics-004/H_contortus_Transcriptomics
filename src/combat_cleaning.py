#!/usr/bin/env python3

"""
Combat-Seq pipeline for raw count matrices.
-------------------------------------------

This script:
- Loads a raw counts matrix
- Automatically generates minimal metadata for batch correction
- Applies Combat-Seq (pycombat_seq)
- Filters corrected counts using CPM > U in at least 3 samples
- Saves the corrected + filtered matrix
- Logs all events and errors
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
    """
    try:
        df = pd.read_csv(path, sep="\t", header=0, index_col=0)
    except Exception as e:
        logging.error(f"Failed to read counts file {path}: {e}")
        raise

    if df.shape[0] < 1 or df.shape[1] < 1:
        logging.error(f"Counts matrix seems empty or malformed: {path} shape={df.shape}")
        raise ValueError("Counts matrix must contain at least 1 gene row and 1 sample column.")

    logging.info(f"Loaded counts matrix '{path}' with shape {df.shape}")
    return df

def parse_batch_flags(batch_flags: List[str], n_columns: int) -> Dict[str, List[int]]:
    """
    Parse repeated --batch flags of format: name:col1,col2,...
    Column numbers are 1-based and include the gene ID column (1 = gene ID).
    Returns mapping {batch_name: [col_numbers (1-based)]}.
    """
    batches = {}
    for entry in batch_flags:
        if ":" not in entry:
            logging.error(f"Invalid --batch entry (missing ':'): {entry}")
            raise ValueError(f"Invalid --batch entry: {entry}. Expected format name:1,2,3")
        name, cols_str = entry.split(":", 1)
        name = name.strip()
        if not name:
            raise ValueError(f"Batch name empty in entry: {entry}")
        cols = [c.strip() for c in cols_str.split(",") if c.strip() != ""]
        if not cols:
            raise ValueError(f"No columns listed for batch '{name}' in entry: {entry}")
        try:
            cols_int = [int(c) for c in cols]
        except ValueError:
            raise ValueError(f"Column indices must be integers in entry: {entry}")
        for ci in cols_int:
            if ci < 1 or ci > n_columns:
                raise ValueError(f"Column index {ci} out of range (1..{n_columns}) for entry: {entry}")
        if name in batches:
            raise ValueError(f"Duplicate batch name: {name}")
        batches[name] = cols_int
    logging.info(f"Parsed batches: {batches}")
    return batches

def build_sample_to_batch_map(batches: Dict[str, List[int]], header_cols: List[str]) -> Dict[str, str]:
    """
    Convert batch column indices (1-based) to sample names and return sample->batch map.

    header_cols: header line split (first element = gene id header, header_cols[1:] = sample names).
    """
    n_columns = len(header_cols)
    sample_names = header_cols[1:]  # columns starting at index 1 correspond to col nums 2..n_columns
    colnum_to_sample = {i + 2: sample_names[i] for i in range(len(sample_names))}  # map 1-based colnums
    sample_to_batch = {}
    for batch_name, colnums in batches.items():
        if len(colnums) == 0:
            raise ValueError(f"Batch '{batch_name}' has no columns assigned.")
        for cn in colnums:
            if cn == 1:
                raise ValueError("Column 1 is the gene ID column; it cannot be assigned to a batch.")
            if cn not in colnum_to_sample:
                raise ValueError(f"Column number {cn} is not a sample column (valid sample columns are 2..{n_columns})")
            sample = colnum_to_sample[cn]
            if sample in sample_to_batch:
                raise ValueError(f"Sample '{sample}' assigned to multiple batches (conflict).")
            sample_to_batch[sample] = batch_name
    logging.info(f"Built sample->batch map for {len(sample_to_batch)} samples")
    return sample_to_batch

def make_metadata_df(sample_to_batch: Dict[str, str], all_samples: List[str]) -> pd.DataFrame:
    """
    Build metadata DataFrame with at least a 'batch' column.
    Unassigned samples get 'unassigned' batch label.
    """
    rows = []
    for s in all_samples:
        batch = sample_to_batch.get(s, "unassigned")
        rows.append((s, batch))
    meta = pd.DataFrame(rows, columns=["sample", "batch"]).set_index("sample")
    logging.info(f"Generated metadata dataframe with shape {meta.shape}")
    return meta

def apply_combat_to_subset(counts: pd.DataFrame, samples_to_correct: List[str], batch_series: pd.Series) -> pd.DataFrame:
    """
    Apply ComBat-seq (negative binomial regression) to a subset of samples.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix (genes x samples).
    samples_to_correct : List[str]
        List of sample column names to apply ComBat-seq to.
    batch_series : pd.Series
        Batch labels for each sample. Index must match columns in counts.
    
    Returns
    -------
    pd.DataFrame
        Corrected counts matrix (same shape as input), with only specified sample columns adjusted.
    """

    for s in samples_to_correct:
        if s not in counts.columns:
            raise ValueError(f"Sample '{s}' not found in counts matrix.")
    
    # Subset counts
    sub_counts = counts.loc[:, samples_to_correct]

     # Align batch vector
    try:
        batches_for_sub = batch_series.loc[samples_to_correct]
    except Exception as e:
        logging.error(f"Failed to align batch vector to selected samples: {e}")
        raise

    try:
        corrected_sub = pycombat_seq(counts=sub_counts, batch=batches_for_sub.tolist())
    except Exception as e:
        logging.error(f"pycombat_seq failed: {e}")
        raise

    corrected_full = counts.copy()
    corrected_full.loc[:, samples_to_correct] = corrected_sub

    logging.info(f"Applied ComBat-seq to {len(samples_to_correct)} samples")
    return corrected_full


def compute_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    libsize = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(libsize, axis=1) * 1e6
    return cpm.fillna(0.0)

def filter_by_cpm(counts: pd.DataFrame, U: float, min_samples: int) -> pd.DataFrame :
    cpm = compute_cpm(counts)
    mask = (cpm > U).sum(axis=1) >= min_samples
    filtered = counts.loc[mask]
    logging.info(f"CPM filter (U={U}, min_samples={min_samples}): before={counts.shape[0]} genes, after={filtered.shape[0]} genes")
    return filtered, cpm

def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    logpath = os.path.join(args.outdir, "combat_seq_pipeline.log")
    configure_logging(logpath)

    if not file_exists_nonempty(args.counts):
        sys.exit("Counts file missing or empty.")
    if not check_matrix_format(args.counts):
        sys.exit("Counts file failed format checks.")

    counts_df = load_counts(args.counts)

    with open(args.counts, "r") as fh:
        header_cols = fh.readline().strip().split("\t")

    batches = parse_batch_flags(args.batch, len(header_cols))
    sample_to_batch = build_sample_to_batch_map(batches, header_cols)
    metadata_df = make_metadata_df(sample_to_batch, list(counts_df.columns))
    metadata_path = os.path.join(args.outdir, "metadata.tsv")
    metadata_df.to_csv(metadata_path, sep="\t")
    logging.info(f"Metadata written to {metadata_path}")

    samples_to_correct = list(sample_to_batch.keys())
    batch_sizes = [len(v) for v in batches.values()]
    min_samples = args.min_samples or min(batch_sizes)
    corrected_df = apply_combat_to_subset(counts_df, samples_to_correct, metadata_df["batch"])

    corrected_path = os.path.join(args.outdir, "corrected_counts.tsv")
    corrected_df.to_csv(corrected_path, sep="\t", float_format="%.6f")
    logging.info(f"Corrected counts saved to {corrected_path}")

    filtered_counts, cpm_matrix = filter_by_cpm(corrected_df, U=args.U, min_samples=min_samples)
    filtered_path = os.path.join(args.outdir, "filtered_counts.tsv")
    filtered_counts.to_csv(filtered_path, sep="\t", float_format="%.6f")
    cpm_path = os.path.join(args.outdir, "cpm_matrix.tsv")
    cpm_matrix.to_csv(cpm_path, sep="\t", float_format="%.6f")
    logging.info(f"Filtered counts saved to {filtered_path}")
    logging.info(f"CPM matrix saved to {cpm_path}")

    print("Done.")
    print(f"Corrected counts: {corrected_path}")
    print(f"Filtered counts:  {filtered_path}")
    print(f"CPM matrix:       {cpm_path}")
    print(f"Metadata:         {metadata_path}")
    print(f"Log:              {logpath}")


if __name__ == "__main__":
    main()
