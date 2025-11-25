#!/usr/bin/env python3
# Standard library imports
import argparse
import logging
import os
import re
import sys
import time
from collections import Counter
from math import isfinite

# Third-party imports
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator
import numpy as np
import pandas as pd
import seaborn as sns
from Bio.KEGG import REST
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Matplotlib configuration
matplotlib.use('Agg')

# ============================
# PARSE CONFIG
# ============================
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Pipeline for KEGG/GO enrichment from DE genes"
    )
    parser.add_argument("--ko_file", required=True, help="TSV with columns: gene_id, KO")
    parser.add_argument("--go_file", required=True, help="TSV with columns: gene_id, GO")
    parser.add_argument("--de_file", required=True, help="DE results CSV with columns: gene_id, log2FoldChange, padj")
    parser.add_argument("--background", required=True, help="Background gene list (one per line)")
    parser.add_argument("--obo_file", required=True, help="GO ontology OBO file")
    parser.add_argument("--outdir", default=None, help="Output directory for all results")
    parser.add_argument("--logfile", default="pipeline.log", help="Log file path")
    parser.add_argument("--fc_thresholds", nargs='+', type=float, default=[1.0], help="log2 fold change thresholds")
    parser.add_argument("--alpha_values", nargs='+', type=float, default=[0.05], help="Adjusted p-value thresholds")
    parser.add_argument("--top_n", type=int, default=20, help="Number of top pathways/GO terms to plot")
    parser.add_argument("--save_files", nargs='+', default=["gene_names","GO_table","KEGG_table","plots"],
                        help="Which files to save: gene_names, GO_table, KEGG_table, plots")
    parser.add_argument("--organism", default="cel", help="KEGG organism code (default: 'cel' for C. elegans)")
    return parser.parse_args()



# ============================
# LOGGING CONFIGURATION
# ============================
def configure_logging(logpath: str):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(logpath), logging.StreamHandler(sys.stdout)]
    )
    logging.info(f"Logging initialized. Logfile: {logpath}")

# ============================
# UTILITIES
# ============================
def file_exists_nonempty(path: str) -> bool:
    """
    Check if file exists and is not empty.
    
    Args:
        path: Path to file to check
        
    Returns:
        True if file exists and has content, False otherwise
    """
    if not os.path.exists(path):
        logging.error(f"File not found: {path}")
        return False
    if os.path.getsize(path) == 0:
        logging.error(f"File is empty: {path}")
        return False
    return True

def validate_inputs(files: list[str]) -> None:
    """
    Validate that all input files exist and are not empty.
    
    Args:
        files: List of file paths to validate
        
    Raises:
        SystemExit: If any file is missing or empty
    """
    for f in files:
        if not file_exists_nonempty(f):
            logging.error(f"Aborting due to missing or empty file: {f}")
            sys.exit(1)
    logging.info("All input files validated.")

def load_inputs(args) -> tuple[pd.DataFrame, set, pd.DataFrame, pd.DataFrame]:
    """
    Load all input files required for the analysis.
    
    Args:
        args: Command line arguments containing file paths
        
    Returns:
        tuple: Contains KO annotations, background genes, GO annotations, and DE results
    """
    logging.info("Loading input files...")
    ko_df = pd.read_csv(args.ko_file, sep="\t", names=['gene_id','KO'])
    background = set(pd.read_csv(args.background, header=None)[0])
    go_df = pd.read_csv(args.go_file, sep="\t", header=None, names=['gene_id','GO'])
    de_df = pd.read_csv(args.de_file, sep="\t", header=0)
    de_df.columns = [c.strip() for c in de_df.columns]
    de_df = de_df.rename(columns={'Unnamed: 0': 'gene_id'})
    print(de_df.columns.tolist())
    return ko_df, background, go_df, de_df


def save_dataframe(df: pd.DataFrame, outdir: str, filename: str, save_files: list[str], key: str) -> None:
    """
    Save DataFrame to file if the key is in save_files list.
    
    Args:
        df: DataFrame to save
        outdir: Output directory path
        filename: Output filename
        save_files: List of file types to save
        key: Key identifying the file type
    """
    if key not in save_files:
        return

    # Create output directory if it doesn't exist
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        path = os.path.join(outdir, filename)
    else:
        path = filename  # Save in current directory if no outdir

    df.to_csv(path, sep="\t", index=False)
    logging.info(f"Saved {key} to {path}")

# ============================
# GENE NAME MAPPING
# ============================
def get_gene_name_from_ko(ko_id: str) -> tuple[str, str | None]:
    """
    Retrieve gene name and symbol from KEGG database using KO identifier.
    
    Args:
        ko_id: KEGG Orthology identifier
        
    Returns:
        tuple: (gene_name, gene_symbol) - Gene name and symbol from KEGG entry
    """
    try:
        entry = REST.kegg_get(ko_id).read()
        name = ko_id
        symbol = None
        
        # Parse KEGG entry to extract name and symbol
        for line in entry.split("\n"):
            if line.startswith("NAME"):
                # Extract gene name from NAME field
                name = line.replace("NAME", "").strip().split(";")[0]
            elif line.startswith("GENES"):
                # Extract gene symbol from GENES field using regex
                # GENES  HSA: 2194(FASN) ...
                match = re.search(r'\(([^)]+)\)', line)
                if match:
                    symbol = match.group(1)

        # Log warnings if information not found
        if name == ko_id:
            logging.warning(f"NAME not found for KO {ko_id}")
        if symbol is None:
            logging.warning(f"GENE SYMBOL not found for KO {ko_id}")

        return name, symbol

    except Exception as e:
        logging.warning(f"Error retrieving KO {ko_id}: {e}")
        return ko_id, None

def map_gene_names_to_new_df(ko_df: pd.DataFrame, outdir: str, save_files: list[str]) -> pd.DataFrame:
    """
    Map KEGG gene names and symbols to KO identifiers for all genes.
    
    Args:
        ko_df: DataFrame with KO identifiers
        outdir: Output directory path
        save_files: List specifying which files to save
        
    Returns:
        DataFrame with gene_id, KO, gene_name, and gene_symbol columns
    """
    logging.info("Mapping gene names from KEGG...")

    gene_names = []
    gene_symbols = []

    # Process each KO identifier to get name and symbol
    for i, ko in enumerate(ko_df.get('KO', [])):
        name, symbol = get_gene_name_from_ko(ko)
        gene_names.append(name)
        gene_symbols.append(symbol)
        
        # Add delay to avoid overwhelming KEGG API
        time.sleep(0.2)

        # Log progress every 50 genes
        if (i + 1) % 50 == 0:
            logging.info(f"Processed {i+1}/{len(ko_df)} genes")

    # Create new DataFrame with mapped names and symbols
    names_df = pd.DataFrame({
        'gene_id': ko_df['gene_id'],
        'KO': ko_df['KO'],
        'gene_name': gene_names,
        'gene_symbol': gene_symbols
    })

    # Save mapping if requested
    save_dataframe(names_df, outdir, "gene_id_to_gene_name.tsv", save_files, "gene_names")
    
    return names_df

# ============================
# DE SPLIT
# ============================
def split_up_down(de_df: pd.DataFrame, fc_thresholds: list[float] = [1.0], 
                  alpha_values: list[float] = [0.05], fc_col: str = 'log2FoldChange', 
                  pval_col: str = 'padj') -> list[dict]:
    """
    Split differentially expressed genes into multiple sets based on FC and p-value thresholds.
    
    Args:
        de_df: DataFrame with DE results
        fc_thresholds: List of fold-change thresholds
        alpha_values: List of p-value thresholds corresponding to fc_thresholds
        fc_col: Column name for fold-change values
        pval_col: Column name for p-values
        
    Returns:
        List of dictionaries containing up/down gene sets for each threshold
    """
    results = []
    # Remove rows with missing FC or p-values
    de_sub = de_df.dropna(subset=[fc_col, pval_col])

    # Ensure thresholds are sorted
    fc_thresholds_sorted = sorted(fc_thresholds)
    
    # Handle single threshold case
    if len(fc_thresholds_sorted) == 1:
        fc = fc_thresholds_sorted[0]
        alpha = alpha_values[0] if alpha_values else 0.05
        
        # Get up and down regulated genes for single threshold
        up = set(de_sub.loc[(de_sub[fc_col] >= fc) & (de_sub[pval_col] < alpha), 'gene_id'])
        down = set(de_sub.loc[(de_sub[fc_col] <= -fc) & (de_sub[pval_col] < alpha), 'gene_id'])
        
        logging.info(f"Single set: FC >= {fc}, alpha={alpha} -> Up: {len(up)}, Down: {len(down)}")
        results.append({'fc': fc, 'alpha': alpha, 'up': up, 'down': down})
    else:
        # Multiple thresholds: create sets for each threshold range
        for i, (fc, alpha) in enumerate(zip(fc_thresholds_sorted, alpha_values)):
            if i == 0:
                # First set: genes between first and second threshold
                upper_fc = fc_thresholds_sorted[i+1] if i+1 < len(fc_thresholds_sorted) else np.inf
                up = set(de_sub.loc[(de_sub[fc_col] >= fc) & (de_sub[fc_col] < upper_fc) & (de_sub[pval_col] < alpha), 'gene_id'])
                down = set(de_sub.loc[(de_sub[fc_col] <= -fc) & (de_sub[fc_col] > -upper_fc) & (de_sub[pval_col] < alpha), 'gene_id'])
            else:
                # Subsequent sets: genes above previous threshold
                lower_fc = fc_thresholds_sorted[i]
                up = set(de_sub.loc[(de_sub[fc_col] >= lower_fc) & (de_sub[pval_col] < alpha), 'gene_id'])
                down = set(de_sub.loc[(de_sub[fc_col] <= -lower_fc) & (de_sub[pval_col] < alpha), 'gene_id'])
            
            logging.info(f"Set {i+1}: FC >= {fc}, alpha={alpha} -> Up: {len(up)}, Down: {len(down)}")
            results.append({'fc': fc, 'alpha': alpha, 'up': up, 'down': down})
    
    return results

def extract_up_down_kos_gos(up_genes: set, down_genes: set, ko_df: pd.DataFrame, 
                           go_df: pd.DataFrame, outdir: str, save_files: list[str]) -> dict:
    """
    Extract KO and GO identifiers for up and down regulated genes.
    
    Args:
        up_genes: Set of up-regulated gene IDs
        down_genes: Set of down-regulated gene IDs
        ko_df: DataFrame with KO annotations
        go_df: DataFrame with GO annotations
        outdir: Output directory path
        save_files: List specifying which files to save
        
    Returns:
        Dictionary with KO and GO sets for up and down regulated genes
    """
    # Extract unique KO identifiers for up and down regulated genes
    up_kos = ko_df.loc[ko_df['gene_id'].isin(up_genes), 'KO'].dropna().unique()
    down_kos = ko_df.loc[ko_df['gene_id'].isin(down_genes), 'KO'].dropna().unique()
    
    # Extract unique GO identifiers for up and down regulated genes
    up_gos = go_df.loc[go_df['gene_id'].isin(up_genes), 'GO'].dropna().unique()
    down_gos = go_df.loc[go_df['gene_id'].isin(down_genes), 'GO'].dropna().unique()

    # Save KO tables if requested
    if "KEGG_table" in save_files:
        save_dataframe(pd.Series(up_kos, name='KO'), outdir, "up_genes_KO.tsv", save_files, "KEGG_table")
        save_dataframe(pd.Series(down_kos, name='KO'), outdir, "down_genes_KO.tsv", save_files, "KEGG_table")
    
    # Save GO tables if requested
    if "GO_table" in save_files:
        save_dataframe(pd.Series(up_gos, name='GO'), outdir, "up_genes_GO.tsv", save_files, "GO_table")
        save_dataframe(pd.Series(down_gos, name='GO'), outdir, "down_genes_GO.tsv", save_files, "GO_table")
    
    return {'up_KO': up_kos, 'down_KO': down_kos, 'up_GO': up_gos, 'down_GO': down_gos}


def go_enrichment(genes_set: set, go_df: pd.DataFrame, de_df: pd.DataFrame, 
                  background: set, obo_file: str, prefix: str, outdir: str = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Perform Gene Ontology enrichment analysis.
    
    Args:
        genes_set: Set of gene IDs for enrichment analysis
        go_df: DataFrame with GO annotations
        de_df: DataFrame with DE results
        background: Set of background gene IDs
        obo_file: Path to GO ontology OBO file
        prefix: Prefix for output files
        outdir: Output directory path
        
    Returns:
        Tuple of (enrichment results DataFrame, extended gene-term mapping DataFrame)
    """
    # Prepare gene to GO mapping (remove NaNs)
    gene2go = go_df.groupby('gene_id')['GO'].apply(lambda s: set([x for x in s if pd.notna(x)])).to_dict()

    # Load GO ontology
    obodag = GODag(obo_file)

    # Convert background to list of strings (GOA requirement)
    population = list(map(str, background))

    # Initialize GO enrichment analysis
    goea = GOEnrichmentStudy(population, gene2go, obodag,
                             propagate_counts=True, alpha=0.05, methods=['fdr_bh'])

    # Run enrichment study (convert genes to strings)
    study_genes = list(map(str, genes_set))
    results_all = goea.run_study(study_genes)

    rows = []
    extended = []

    # Process significant GO terms
    for r in results_all:
        # Handle different attribute names for FDR across versions
        p_fdr = getattr(r, "p_fdr_bh", None)
        if p_fdr is None:
            p_fdr = getattr(r, "p_fdr", None) or getattr(r, "pvalue_fdr_bh", None)

        if p_fdr is None:
            continue

        # Only consider significant terms (FDR < 0.05)
        if p_fdr < 0.05:
            # Get namespace safely across different versions
            ns = getattr(r, "NS", None)
            if ns is None:
                ns = getattr(r, "namespace", None)
            # Handle different namespace formats
            if isinstance(ns, (dict, list)):
                try:
                    if isinstance(ns, dict):
                        ns = ns.get("namespace") or ns.get("id") or ns.get("name")
                    else:
                        ns = str(ns)
                except Exception:
                    ns = str(ns)
            # Map full names to standard codes (BP/MF/CC)
            ns_map_full = {
                'biological_process': 'BP',
                'molecular_function': 'MF',
                'cellular_component': 'CC',
                'Biological Process': 'BP',
                'Molecular Function': 'MF',
                'Cellular Component': 'CC'
            }
            ns_code = ns_map_full.get(str(ns), str(ns))

            # Calculate enrichment ratios safely
            ratio_in_study = np.nan
            ratio_in_pop = np.nan
            try:
                if getattr(r, 'ratio_in_study', None):
                    ratio_in_study = r.ratio_in_study[0] / r.ratio_in_study[1]
                if getattr(r, 'ratio_in_pop', None):
                    ratio_in_pop = r.ratio_in_pop[0] / r.ratio_in_pop[1]
            except Exception:
                pass

            # Add to results table
            rows.append({
                'GO': r.GO,
                'name': getattr(r, 'name', ''),
                'namespace': ns_code,
                'p_fdr': float(p_fdr),
                'ratio_in_study': ratio_in_study,
                'ratio_in_pop': ratio_in_pop,
                'n_genes': len(r.study_items) if getattr(r, 'study_items', None) else 0
            })

            # Extended table: list genes in study for each term
            study_items = getattr(r, 'study_items', [])
            genes_in_term = [g for g in study_items if str(g) in set(map(str, genes_set))]
            for g in genes_in_term:
                # Look up gene in DE results if available
                row_df = None
                try:
                    mask = de_df['gene_id'].astype(str) == str(g)
                    if mask.any():
                        row_df = de_df.loc[mask].iloc[0]
                except Exception:
                    row_df = None
                extended.append({
                    'GO': r.GO,
                    'GO_name': getattr(r, 'name', ''),
                    'namespace': ns_code,
                    'Gene_ID': str(g),
                    'Gene_Name': row_df.get('gene_name', str(g)) if row_df is not None and 'gene_name' in row_df else str(g),
                    'log2FC': row_df.get('log2FoldChange', np.nan) if row_df is not None and 'log2FoldChange' in row_df else np.nan,
                    'p_fdr': float(p_fdr),
                    'pvalue': row_df.get('padj', np.nan) if row_df is not None and 'padj' in row_df else np.nan,
                    'ratio_in_study': ratio_in_study,
                    'ratio_in_pop': ratio_in_pop
                })

    results_df = pd.DataFrame(rows)
    extended_df = pd.DataFrame(extended)

    # Save output files
    fname1 = f"{prefix}_GO_enrichment.tsv"
    fname2 = f"{prefix}_GO_extended_table.tsv"
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        results_df.to_csv(os.path.join(outdir, fname1), sep="\t", index=False)
        extended_df.to_csv(os.path.join(outdir, fname2), sep="\t", index=False)
    else:
        results_df.to_csv(fname1, sep="\t", index=False)
        extended_df.to_csv(fname2, sep="\t", index=False)

    logging.info(f"{len(results_df)} significant GO terms saved for {prefix}")
    logging.info(f"Extended GO table with genes saved: {len(extended_df)} rows")

    # Generate plots if results exist
    if not results_df.empty:
        plot_go_enrichment(results_df, prefix, top_n=20, outdir=outdir)
    else:
        logging.warning(f"No significant GO terms found for {prefix}")

    return results_df, extended_df

def plot_go_enrichment(results: pd.DataFrame, prefix: str, top_n: int = 20, outdir: str = None):
    """
    Generate GO enrichment visualization plots.
    
    Args:
        results: DataFrame with GO enrichment results
        prefix: Prefix for output files
        top_n: Number of top terms to plot
        outdir: Output directory path
    """
    if results is None or results.empty:
        print(f"No results for {prefix}")
        return

    # Normalize namespace to BP/MF/CC codes
    def _norm_ns(x):
        if pd.isna(x):
            return None
        if isinstance(x, (dict, list)):
            try:
                if isinstance(x, dict):
                    return x.get('namespace') or x.get('id') or str(x)
                return str(x)
            except Exception:
                return str(x)
        return str(x)

    results = results.copy()
    results['namespace'] = results['namespace'].apply(_norm_ns).astype(str)
    ns_map = {
        'biological_process': 'BP',
        'molecular_function': 'MF',
        'cellular_component': 'CC',
        'BP': 'BP', 'MF': 'MF', 'CC': 'CC'
    }
    results['namespace'] = results['namespace'].map(lambda x: ns_map.get(x, x))

    # Select top N terms by FDR
    top = results.sort_values('p_fdr').head(top_n).copy()
    
    # Debug information
    print(f"\n=== GO RESULTS {prefix} ===")
    print(f"Total terms: {len(results)}, Top terms: {len(top)}")
    print(f"FDR range: {results['p_fdr'].min():.2e} to {results['p_fdr'].max():.2e}")
    print(f"Ratio_in_study range: {results['ratio_in_study'].min():.4f} to {results['ratio_in_study'].max():.4f}")
    print(f"Ratio_in_pop range: {results['ratio_in_pop'].min():.4f} to {results['ratio_in_pop'].max():.4f}")
    
    # Calculate enrichment ratio
    top['enrichment_ratio'] = top['ratio_in_study'] / top['ratio_in_pop']
    top['enrichment_ratio'] = top['enrichment_ratio'].replace([np.inf, -np.inf], np.nan)
    
    print(f"Enrichment ratio range: {top['enrichment_ratio'].min():.4f} to {top['enrichment_ratio'].max():.4f}")
    print(f"Gene counts range: {results['n_genes'].min()} to {results['n_genes'].max()}")
    
    # Show top 5 terms with their ratios
    print("\nTop 5 terms with ratios:")
    for i, row in top.head().iterrows():
        print(f"  {row['name'][:50]}... | FDR: {row['p_fdr']:.2e} | Study: {row['ratio_in_study']:.4f} | Pop: {row['ratio_in_pop']:.4f} | Enrich: {row.get('enrichment_ratio', 'NaN'):.4f}")
    
    # Calculate -log10(FDR) for visualization
    top['neg_log10_fdr'] = -np.log10(top['p_fdr'].replace(0, np.nextafter(0, 1)))

    # Assign colors by namespace
    ns_colors = {'BP': "#1d96ed", 'MF': "#f62020", 'CC': "#24e824"}
    top['color'] = top['namespace'].map(lambda x: ns_colors.get(x, '#7f7f7f'))

    # Truncate long names for display
    top['short_name'] = top['name'].apply(lambda s: (s[:80] + '...') if isinstance(s, str) and len(s) > 80 else s)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

    # --- PLOT 1: Horizontal bars (-log10 FDR) ---
    if len(top) > 0:
        bars = ax1.barh(y=np.arange(len(top)), width=top['neg_log10_fdr'], color=top['color'])
        ax1.set_yticks(np.arange(len(top)))
        ax1.set_yticklabels(top['short_name'], fontsize=9)
        ax1.set_xlabel('-log10(FDR)')
        ax1.set_title(f"Top {len(top)} GO Terms — {prefix}", fontweight='bold')
        
        # Adjust X-axis limits for bar plot
        max_fdr = top['neg_log10_fdr'].max()
        ax1.set_xlim(0, max_fdr * 1.15)
        
        print(f"Bar plot X limits: 0 to {max_fdr * 1.15:.2f}")
        
        # Add annotations with gene count and namespace
        for i, (_, row) in enumerate(top.iterrows()):
            txt = f"n={int(row.get('n_genes', 0))} | {row.get('namespace','')}"
            ax1.text(row['neg_log10_fdr'] + max_fdr * 0.01, i, txt, va='center', fontsize=8)
    else:
        ax1.text(0.5, 0.5, 'No significant terms', 
                transform=ax1.transAxes, ha='center', va='center', fontsize=12)
        ax1.set_title(f"GO Terms — {prefix}")

    # --- PLOT 2: Bubble plot (Enrichment Ratio) ---
    # Filter out NaN values and reset indices for bubble plot
    bubble_data = top.dropna(subset=['enrichment_ratio']).copy().reset_index(drop=True)
    
    if len(bubble_data) > 0:
        # Calculate bubble sizes based on gene counts
        bubble_sizes = (bubble_data['n_genes'].fillna(1).astype(float)) * 25
        
        sc = ax2.scatter(bubble_data['enrichment_ratio'], np.arange(len(bubble_data)), 
                        s=bubble_sizes,
                        c=bubble_data['neg_log10_fdr'],
                        cmap='viridis', alpha=0.8, edgecolor='k')
        
        ax2.set_yticks(np.arange(len(bubble_data)))
        ax2.set_yticklabels(bubble_data['short_name'], fontsize=9)
        
        # Adjust X-axis limits for bubble plot
        min_ratio = bubble_data['enrichment_ratio'].min()
        max_ratio = bubble_data['enrichment_ratio'].max()
        
        # Ensure adequate margins
        x_margin = (max_ratio - min_ratio) * 0.1
        if x_margin == 0:  # If all values are equal
            x_margin = max_ratio * 0.1 if max_ratio != 0 else 0.1
            
        ax2.set_xlim(min_ratio - x_margin, max_ratio + x_margin)
        
        # Debug: show applied limits
        print(f"Bubble plot X limits: {min_ratio - x_margin:.3f} to {max_ratio + x_margin:.3f}")
        
        # Add colorbar
        cbar = plt.colorbar(sc, ax=ax2)
        cbar.set_label('-log10(FDR)')
        
        # Reference line at enrichment ratio = 1
        ax2.axvline(1.0, color='red', linestyle='--', alpha=0.6)
    else:
        ax2.text(0.5, 0.5, 'No enrichment ratio data', 
                transform=ax2.transAxes, ha='center', va='center', fontsize=12)
        print("No enrichment ratio data for bubble plot")
    
    ax2.invert_yaxis()
    ax2.set_xlabel('Enrichment Ratio (Study / Population)')
    ax2.set_title('GO Enrichment Bubble Plot')

    plt.tight_layout()
    outpath = f"{prefix}_GO_enrichment.png"
    if outdir:
        outpath = os.path.join(outdir, outpath)
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")

    # Generate small panel by namespace
    plot_go_by_namespace(results, prefix, top_n=10, outdir=outdir)

def plot_go_by_namespace(results: pd.DataFrame, prefix: str, top_n: int = 10, outdir: str = None):
    """
    Generate GO enrichment plots separated by namespace (BP, MF, CC).
    
    Args:
        results: DataFrame with GO enrichment results
        prefix: Prefix for output files
        top_n: Number of top terms to plot per namespace
        outdir: Output directory path
    """
    if results is None or results.empty:
        return

    # Normalize namespace (same logic as before)
    def _norm_ns(x):
        if pd.isna(x):
            return None
        if isinstance(x, (dict, list)):
            try:
                if isinstance(x, dict):
                    return x.get('namespace') or x.get('id') or str(x)
                return str(x)
            except Exception:
                return str(x)
        return str(x)

    results = results.copy()
    results['namespace'] = results['namespace'].apply(_norm_ns).astype(str)
    ns_map = {
        'biological_process': 'BP',
        'molecular_function': 'MF',
        'cellular_component': 'CC',
        'BP': 'BP', 'MF': 'MF', 'CC': 'CC'
    }
    results['namespace'] = results['namespace'].map(lambda x: ns_map.get(x, x))

    namespaces = ['BP', 'MF', 'CC']
    colors = {'BP': '#1f77b4', 'MF': '#ff7f0e', 'CC': '#2ca02c'}

    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(22, 6), sharey=False)

    for ax, ns in zip(axes, namespaces):
        # Get top terms for current namespace
        sub = results[results['namespace'] == ns].sort_values('p_fdr').head(top_n).copy()
        if sub is None or sub.empty:
            ax.set_title(f"{ns}: no terms")
            ax.axis('off')
            continue
        
        # Calculate -log10(FDR) and truncate names
        sub['neg_log10_fdr'] = -np.log10(sub['p_fdr'].replace(0, np.nextafter(0, 1)))
        sub['short_name'] = sub['name'].apply(lambda s: (s[:50] + '...') if isinstance(s, str) and len(s) > 50 else s)
        
        # Create horizontal bar plot
        bars = ax.barh(y=np.arange(len(sub)), width=sub['neg_log10_fdr'], color=colors.get(ns, '#7f7f7f'))
        
        # Adjust X-axis limits for each subplot individually
        max_fdr_sub = sub['neg_log10_fdr'].max()
        ax.set_xlim(0, max_fdr_sub * 1.15)
        
        # Set Y-axis labels and formatting
        ax.set_yticks(np.arange(len(sub)))
        ax.set_yticklabels(sub['short_name'], fontsize=8)
        ax.set_xlabel('-log10(FDR)')
        ax.set_title(f"{ns} (n={len(sub)})")
        ax.invert_yaxis()

    plt.tight_layout()
    outpath = f"{prefix}_GO_by_namespace_corrected.png"
    if outdir:
        outpath = os.path.join(outdir, outpath)
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()


def _chunked_iterable(iterable, size: int):
    """
    Split an iterable into chunks of specified size.
    
    Args:
        iterable: Input iterable to chunk
        size: Size of each chunk
        
    Yields:
        List of chunks from the original iterable
    """
    it = list(iterable)
    # Yield chunks of specified size
    for i in range(0, len(it), size):
        yield it[i:i+size]

def _parse_kegg_entry_text(entry_text: str) -> tuple[str | None, list[str], str]:
    """
    Parse a single KEGG entry block text and extract pathway information.
    
    Args:
        entry_text: Raw text from KEGG REST response
        
    Returns:
        Tuple of (entry_id, pathway_list, title)
    """
    entry_id = None
    paths = []
    current_section = None
    last_path = None
    title = ""
    
    # Parse each line of the KEGG entry
    for line in entry_text.splitlines():
        if line.startswith("ENTRY"):
            # Extract entry ID from ENTRY line
            parts = line.split()
            if len(parts) >= 2:
                entry_id = parts[1]
        elif line.startswith("NAME"):
            # Extract entry title from NAME line
            title = line.replace("NAME", "").strip()
        elif line.startswith("PATHWAY"):
            # Extract pathway IDs from PATHWAY lines
            # PATHWAY   map00010  Glycolysis / Gluconeogenesis
            parts = line.split()
            if len(parts) >= 2:
                pid = parts[1]
                paths.append(pid)
                last_path = pid
        elif line.startswith("CLASS"):
            # Mark start of CLASS section
            current_section = "CLASS"
        elif current_section == "CLASS" and line.startswith(" "):
            # CLASS section lines are indented
            # Class information is handled externally
            pass
        else:
            # Reset section on blank lines
            if line.strip() == "":
                current_section = None
                
    return entry_id, paths, title

def get_ko_pathway_mapping(ko_list: list[str], chunk_size: int = 20, sleep: float = 0.2, 
                          organism: str = None) -> tuple[dict, dict]:
    """
    Map KO IDs to pathways using KEGG REST efficient multi-fetch.
    
    Args:
        ko_list: List of KO identifiers
        chunk_size: Number of KOs per REST API call
        sleep: Pause between chunk requests in seconds
        organism: KEGG organism code to filter pathways
        
    Returns:
        Tuple of (KO to pathway mapping, pathway information dictionary)
    """
    # Normalize KO IDs as strings and remove duplicates
    kos = [str(k).strip() for k in ko_list if pd.notna(k)]
    kos = list(dict.fromkeys(kos))

    # Normalize KO formats (remove 'ko:' or 'K:' prefixes)
    def _norm_ko(x):
        return re.sub(r'^(ko:|K:)?', '', x, flags=re.IGNORECASE).strip()

    kos_norm = list(map(_norm_ko, kos))

    ko2path = {}
    pathway_info = {}

    # Precompute valid pathway numeric suffixes for organism filtering
    valid_path_digits = None
    if organism:
        try:
            listing = REST.kegg_list("pathway", organism).read()
            valid_path_digits = set()
            for line in listing.splitlines():
                if not line.strip(): 
                    continue
                left = line.split("\t")[0].strip()  # like 'path:cel00010'
                m = re.search(r'(\d+)$', left)
                if m:
                    valid_path_digits.add(m.group(1))
        except Exception as e:
            logging.warning(f"Could not fetch pathway list for organism '{organism}': {e}")
            valid_path_digits = None

    cache = {}  # In-run cache: ko -> paths

    # Process KOs in chunks to avoid overwhelming API
    for chunk in _chunked_iterable(kos_norm, chunk_size):
        query = "+".join(chunk)
        try:
            text = REST.kegg_get(query).read()
        except Exception as e:
            logging.warning(f"KEGG REST multi-get failed for chunk starting with {chunk[0]}: {e}")
            # Fallback to single-get per KO
            text = ""
            for k in chunk:
                try:
                    text += REST.kegg_get(k).read() + "\n\n"
                except Exception as e2:
                    logging.warning(f"KEGG REST single-get failed for {k}: {e2}")
                time.sleep(sleep)
                
        # Split response into individual entry blocks
        blocks = re.split(r'\n(?=ENTRY\s)', text)
        for block in blocks:
            if not block.strip():
                continue
                
            # Parse entry block to extract basic information
            eid, paths, title = _parse_kegg_entry_text(block)
            if not eid:
                # Try to infer entry ID from first line
                firstline = block.splitlines()[0]
                eid = firstline.split()[1] if len(firstline.split())>1 else None
                
            # Choose canonical KO key from the block
            ko_key_candidates = []
            if eid:
                ko_key_candidates.append(re.sub(r'^(ko:|K:)?', '', eid, flags=re.IGNORECASE))
                
            if ko_key_candidates:
                ko_key = ko_key_candidates[0]
            else:
                # Skip block if no valid KO key found
                continue

            # Parse PATHWAY lines more robustly from the block
            block_paths = []
            for line in block.splitlines():
                if line.startswith("PATHWAY"):
                    parts = line.split()
                    if len(parts) >= 2:
                        pid = parts[1]
                        block_paths.append(pid)
                        # Store pathway name if not already present
                        if pid not in pathway_info:
                            pname = " ".join(parts[2:]) if len(parts) > 2 else "Unknown"
                            pathway_info[pid] = {'name': pname, 'class': 'Unknown'}
                            
            # Parse CLASS lines for pathway classification
            cur_class = None
            for i, line in enumerate(block.splitlines()):
                if line.startswith("CLASS"):
                    cur_class = ""
                    # Collect all indented lines in CLASS section
                    j = i+1
                    while j < len(block.splitlines()) and block.splitlines()[j].startswith(" "):
                        cur_class += block.splitlines()[j].strip() + " "
                        j += 1
                    cur_class = cur_class.strip()
                    # Assign class to all pathways in this block
                    for p in block_paths:
                        if p in pathway_info:
                            pathway_info[p]['class'] = cur_class

            # Store pathways for this KO (filtering happens later)
            cache[ko_key] = block_paths

        time.sleep(sleep)

    # Build final ko2path with optional organism filtering
    for k, paths in cache.items():
        if not paths:
            continue
        kept = []
        if valid_path_digits is not None:
            # Filter pathways by organism-specific numeric suffixes
            for p in paths:
                m = re.search(r'(\d+)$', str(p))
                if m and m.group(1) in valid_path_digits:
                    kept.append(p)
        else:
            kept = paths[:]
        if kept:
            ko2path[k] = kept

    # Filter pathway_info to only include pathways present in ko2path
    if pathway_info:
        allowed = set([p for paths in ko2path.values() for p in paths])
        pathway_info = {p: pathway_info[p] for p in pathway_info if p in allowed}

    return ko2path, pathway_info

def fisher_enrichment(ko2path_study: dict, background_genes: list[str], pathway_info: dict, 
                     organism: str = None, chunk_size: int = 20, sleep: float = 0.2) -> pd.DataFrame:
    """
    Compute Fisher enrichment for KEGG pathways.
    
    Args:
        ko2path_study: Dict of study KO to pathway mapping
        background_genes: List of background KOs
        pathway_info: Pathway metadata dictionary
        organism: KEGG organism code for filtering
        chunk_size: Chunk size for API calls
        sleep: Sleep time between API calls
        
    Returns:
        DataFrame with enrichment results and FDR-corrected p-values
    """
    # Normalize background KOs and remove duplicates
    bg_kos = [str(x).strip() for x in background_genes if pd.notna(x)]
    bg_kos_norm = list(dict.fromkeys([re.sub(r'^(ko:|K:)?', '', x, flags=re.IGNORECASE) for x in bg_kos]))

    # Fetch mapping for background KOs not already in study mapping
    missing_bg = [k for k in bg_kos_norm if k not in ko2path_study]
    bg_ko2path = {}
    if missing_bg:
        bg_ko2path, _ = get_ko_pathway_mapping(missing_bg, chunk_size=chunk_size, sleep=sleep, organism=organism)
        
    # Create combined background mapping using study data when available
    combined_bg_map = {}
    for k in bg_kos_norm:
        if k in ko2path_study:
            combined_bg_map[k] = ko2path_study[k]
        elif k in bg_ko2path:
            combined_bg_map[k] = bg_ko2path[k]
        else:
            combined_bg_map[k] = []

    # Count pathways in study and background
    study_counts = Counter([p for paths in ko2path_study.values() for p in paths])
    bg_counts = Counter([p for paths in combined_bg_map.values() for p in paths])

    results = []
    total_study = len(ko2path_study)
    total_bg = len([k for k in combined_bg_map.keys()])

    # Perform Fisher's exact test for each pathway
    for pid, count in study_counts.items():
        bg_count = bg_counts.get(pid, 0)
        # Create contingency table: study vs background in pathway
        table = [[count, total_study - count],
                 [bg_count, total_bg - bg_count]]
        # Calculate Fisher's exact test
        try:
            odds, pval = fisher_exact(table, alternative='greater')
        except Exception as e:
            odds, pval = (np.nan, 1.0)
            
        results.append({
            'Pathway': pid,
            'Pathway_name': pathway_info.get(pid, {}).get('name', 'Unknown'),
            'Pathway_class': pathway_info.get(pid, {}).get('class', 'Unknown'),
            'Study_count': int(count),
            'Background_count': int(bg_count),
            'Odds_ratio': odds,
            'P_value': float(pval)
        })

    df = pd.DataFrame(results)
    # Apply FDR correction if results exist
    if not df.empty:
        df['FDR'] = multipletests(df['P_value'], method='fdr_bh')[1]
        df.sort_values('FDR', inplace=True)
    return df

def plot_kegg_results(df: pd.DataFrame, prefix: str, top_n: int = 20, outdir: str = None):
    """
    Generate KEGG pathway enrichment visualization plots.
    
    Args:
        df: DataFrame with KEGG enrichment results
        prefix: Prefix for output files
        top_n: Number of top pathways to plot
        outdir: Output directory path
    """
    if df is None or df.empty:
        logging.warning(f"No KEGG results to plot for {prefix}")
        return

    df = df.copy()
    # Ensure numeric data types
    df['Study_count'] = pd.to_numeric(df['Study_count'], errors='coerce').fillna(0).astype(int)
    df['FDR'] = pd.to_numeric(df.get('FDR', df.get('P_value')), errors='coerce').fillna(1.0)
    df['neg_log10_fdr'] = -np.log10(df['FDR'].replace(0, np.nextafter(0,1)))

    # Select top pathways and prepare for plotting
    top = df.head(top_n).copy()
    top['short_name'] = top['Pathway_name'].apply(lambda s: (s[:80] + '...') if isinstance(s, str) and len(s) > 80 else s)
    top = top.sort_values('Study_count', ascending=True)  # Largest on top after invert

    # -------------------
    # BAR PLOT
    # -------------------
    fig, ax = plt.subplots(figsize=(10, max(6, 0.3 * len(top))))
    cmap = plt.get_cmap('viridis')
    colors = cmap((top['neg_log10_fdr'] - top['neg_log10_fdr'].min()) / 
                 (top['neg_log10_fdr'].max() - top['neg_log10_fdr'].min() + 1e-9))

    # Create horizontal bar plot
    ax.barh(range(len(top)), top['Study_count'], color=colors)

    # Set Y-axis labels with odds ratio information
    new_labels = [f"{row['short_name']}  (OR={row['Odds_ratio']:.2g})" for _, row in top.iterrows()]
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(new_labels, fontsize=9)

    ax.set_xlabel('Number of KOs (Study_count)')
    ax.set_title(f"Top {len(top)} KEGG Pathways ({prefix})")

    # Add colorbar for significance
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=top['neg_log10_fdr'].min(), 
                                                                   vmax=top['neg_log10_fdr'].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('-log10(FDR)')

    plt.tight_layout()
    outpath = f"{prefix}_KEGG_barplot.png"
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, outpath)
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()

    # -------------------
    # BUBBLE PLOT
    # -------------------
    bubble = top.copy().reset_index(drop=True)
    bubble_fig, bubble_ax = plt.subplots(figsize=(10, max(6, 0.3 * len(bubble))))

    # Scale bubble sizes based on study count
    min_size = 50
    max_size = 1000
    counts = bubble['Study_count'].astype(float)
    sizes = min_size + (counts - counts.min()) / (counts.max() - counts.min() + 1e-9) * (max_size - min_size)

    # Create scatter plot with bubbles
    sc = bubble_ax.scatter(
        bubble['Odds_ratio'].replace([np.inf, -np.inf], np.nan),
        bubble.index,
        s=sizes,
        c=bubble['neg_log10_fdr'],
        cmap='viridis',
        alpha=0.8,
        edgecolor='k'
    )

    # Set Y-axis labels with pathway IDs
    new_labels = [f"{row['short_name']}  ({str(row['Pathway'])[:8]})" for _, row in bubble.iterrows()]
    bubble_ax.yaxis.set_major_locator(FixedLocator(bubble.index))
    bubble_ax.set_yticklabels(new_labels, fontsize=9)

    bubble_ax.set_xlabel('Odds ratio')
    bubble_ax.set_title(f"KEGG Bubble Plot ({prefix})")
    # Add reference line at odds ratio = 1
    bubble_ax.axvline(1.0, color='red', linestyle='--', alpha=0.6)

    # Add colorbar for significance
    cbar2 = plt.colorbar(sc, ax=bubble_ax)
    cbar2.set_label('-log10(FDR)')

    outbubble = f"{prefix}_KEGG_bubbleplot.png"
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        outbubble = os.path.join(outdir, outbubble)

    plt.tight_layout()
    plt.savefig(outbubble, dpi=300, bbox_inches='tight')
    plt.close()

def kegg_enrichment(ko_list: list[str], background_genes: list[str], prefix: str = "UP", 
                   top_n: int = 20, organism: str = "cel", chunk_size: int = 20, 
                   sleep: float = 0.2, outdir: str = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Full KEGG enrichment pipeline with optimized batch processing.
    
    Args:
        ko_list: List of study KOs
        background_genes: List of background KOs
        prefix: Prefix for output files
        top_n: Number of top pathways to plot
        organism: KEGG organism code for pathway filtering
        chunk_size: Chunk size for API calls
        sleep: Sleep time between API calls
        outdir: Output directory path
        
    Returns:
        Tuple of (summary DataFrame, extended mapping DataFrame)
    """
    # Normalize and deduplicate KO identifiers
    study_kos = [re.sub(r'^(ko:|K:)?', '', str(x).strip(), flags=re.IGNORECASE) for x in ko_list if pd.notna(x)]
    bg_kos = [re.sub(r'^(ko:|K:)?', '', str(x).strip(), flags=re.IGNORECASE) for x in background_genes if pd.notna(x)]
    union_kos = list(dict.fromkeys(study_kos + bg_kos))

    if len(union_kos) == 0:
        logging.warning("No KOs provided for KEGG enrichment.")
        return None, None

    # Fetch KO to pathway mapping for all KOs in one pass
    ko2path_all, pathway_info_all = get_ko_pathway_mapping(union_kos, chunk_size=chunk_size, sleep=sleep, organism=organism)

    # Extract study-specific mapping
    ko2path_study = {k: v for k, v in ko2path_all.items() if k in study_kos}
    # Compute Fisher enrichment
    df_summary = fisher_enrichment(ko2path_study, bg_kos, pathway_info_all, organism=organism, chunk_size=chunk_size, sleep=sleep)

    if df_summary is None or df_summary.empty:
        logging.warning(f"No KEGG enrichment results for {prefix}")
        return None, None

    # Save summary results
    fname_summary = f"{prefix}_KEGG_summary.tsv"
    fname_extended = f"{prefix}_KEGG_extended_table.tsv"
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        df_summary.to_csv(os.path.join(outdir, fname_summary), sep="\t", index=False)
    else:
        df_summary.to_csv(fname_summary, sep="\t", index=False)

    # Create extended table with KO-pathway mappings
    extended_rows = []
    for ko, paths in ko2path_study.items():
        for p in paths:
            extended_rows.append({
                'KO': ko,
                'Pathway': p,
                'Pathway_name': pathway_info_all.get(p, {}).get('name', 'Unknown'),
                'Pathway_class': pathway_info_all.get(p, {}).get('class', 'Unknown')
            })
    df_extended = pd.DataFrame(extended_rows)
    if outdir:
        df_extended.to_csv(os.path.join(outdir, fname_extended), sep="\t", index=False)
    else:
        df_extended.to_csv(fname_extended, sep="\t", index=False)

    # Generate visualization plots
    plot_kegg_results(df_summary, prefix, top_n=top_n, outdir=outdir)

    return df_summary, df_extended

# ============================
# MAIN FUNCTION
# ============================

def main():
    """
    Execute KEGG and GO enrichment analysis pipeline.
    """
    # Parse command line arguments
    args = parse_arguments()

    # Configure logging
    configure_logging(args.logfile)

    # Validate input files exist and are not empty
    validate_inputs([args.ko_file, args.go_file, args.de_file, args.background, args.obo_file])

    # Load input data
    ko_df, background, go_df, de_df = load_inputs(args)

    # Optional: Map gene names from KEGG (commented out for performance)
    # names_df = map_gene_names_to_new_df(ko_df, outdir=args.outdir, save_files=args.save_files)

    # Split DE genes into sets based on fold-change and p-value thresholds
    de_sets = split_up_down(
        de_df,
        fc_thresholds=args.fc_thresholds,
        alpha_values=args.alpha_values,
        fc_col='log2FoldChange',
        pval_col='padj'
    )

    # Extract background KOs for enrichment analysis
    background_kos = (
        ko_df[ko_df["gene_id"].isin(background)]["KO"]
        .dropna()
        .unique()
        .tolist()
    )

    if len(background_kos) == 0:
        logging.error("No background genes have KO annotations — cannot run KEGG.")
        sys.exit(1)

    # Process each DE gene set
    for i, s in enumerate(de_sets):
        up_genes = s['up']
        down_genes = s['down']
        
        # Skip empty sets
        if len(up_genes) == 0 and len(down_genes) == 0:
            logging.info(f"Skipping Set{i+1} - no genes found")
            continue

        prefix = f"Set{i+1}_FC{str(s['fc']).replace('.', '_')}_alpha{str(s['alpha']).replace('.', '_')}"
        logging.info(f"Processing {prefix} ...")

        print("UP:", len(up_genes))
        print("DOWN:", len(down_genes))
        print("Background:", len(background))

        # Extract KO and GO identifiers for up/down regulated genes
        files = extract_up_down_kos_gos(up_genes, down_genes, ko_df, go_df, outdir=args.outdir, save_files=args.save_files)

        # Perform GO enrichment for up-regulated genes if present
        if len(up_genes) > 0:
            go_up, go_up_ext = go_enrichment(up_genes, go_df, de_df, background, args.obo_file, prefix=prefix + "_UP", outdir=args.outdir)
        else:
            logging.info(f"No UP genes for GO enrichment in {prefix}")
            
        # Perform GO enrichment for down-regulated genes if present
        if len(down_genes) > 0:
            go_down, go_down_ext = go_enrichment(down_genes, go_df, de_df, background, args.obo_file, prefix=prefix + "_DOWN", outdir=args.outdir)
        else:
            logging.info(f"No DOWN genes for GO enrichment in {prefix}")

        # Perform KEGG enrichment for up-regulated KOs if present
        if len(files['up_KO']) > 0:
            kegg_up, kegg_up_ext = kegg_enrichment(
                ko_list=list(files['up_KO']),
                background_genes=background_kos,
                prefix=prefix + "_UP",
                top_n=args.top_n,
                outdir=args.outdir,
                organism=args.organism
            )
        else:
            logging.info(f"No UP KOs for KEGG enrichment in {prefix}")

        # Perform KEGG enrichment for down-regulated KOs if present
        if len(files['down_KO']) > 0:
            kegg_down, kegg_down_ext = kegg_enrichment(
                ko_list=list(files['down_KO']),
                background_genes=background_kos,
                prefix=prefix + "_DOWN", 
                top_n=args.top_n,
                outdir=args.outdir,
                organism=args.organism
            )
        else:
            logging.info(f"No DOWN KOs for KEGG enrichment in {prefix}")

    logging.info("Pipeline finished successfully.")

if __name__ == "__main__":
    main()
