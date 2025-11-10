#!/usr/bin/env python3
'''
Data visualization before analysis
----------------------------------
Given a count matrix from gene expresion data, allows the visualization of the data before
doing a Differential Gene Expression analysis. It supose the next structure of the count matrix:

gene_id,gene_name,sample1,sample2,...,samplen

Authors: Yael Montiel-Vargas & Pablo Salazar-Méndez
Version: 0.0.0
Date: 2025-11-07

This script seeks to ease the decision making process for a succesful DGE analysis by helping in
the visualizing process with the following pipeline:

    - Restructure the count matrix in function of the conditions that shall be compared.
    - Logarithmic (base 2) normalization of the data.
    - Count values melting within a single column.
    - Data distribution visualization:
        - KDE
        - Dimensional reduction (PCA, MDS or both)
    - Save images in the given direction (default is the current working directory).

Dependencies:
    - matplotlib.pyplot
    - seaborn
    - pandas
    - scikitlearn
    - numpy

Example usage:
    python3 pre_de_visualization.py \
        -m results/count_mtx.csv -c Treated,Untreated \
        --blocks 3,4 -s "," -l 2 -o 4-6,1-3,7 \
        --start-at 1 -O results/Images
'''

from scipy.spatial.distance import pdist, squareform
from argparse import Namespace, ArgumentParser
from pandas import DataFrame, read_csv
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from scipy.stats import zscore
from datetime import datetime

import matplotlib.pyplot as plt
import multiprocessing as mp
import logging as lg
import seaborn as sns
import numpy as np
import sys
import re
import os

# ============================================================================= #
#                                   UTILITY
# ============================================================================= #

log_filename = f"pre_visualization_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
lg.basicConfig(
    level=lg.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[lg.FileHandler(log_filename), lg.StreamHandler()]
)

def my_parser() -> Namespace:
    parser = ArgumentParser(description="Data visualization before DGE analysis pipeline")
    parser.add_argument("-m", "--matrix", required=True, type=str, help="Count matrix direction")
    parser.add_argument("-c", "--conditions", required=True, type=str, help="Condition to be compared separated by a coma (,)")
    parser.add_argument("-b", "--blocks", required=True, type=str, help="Number of blocks associated to each condition")
    parser.add_argument("-s", "--separator", type=str, default="$'\t'", help="Separator of the count matrix (default is a tab)")
    parser.add_argument("-l", "--log-base", type=int, default=2, help="Base for the logarithmic normalization")
    parser.add_argument("-o", "--order-col", required=True, type=str, help="Final order of the columns within the DataFrame (Pandas)")
    parser.add_argument("--start-at", type=int, default=1, help="Sufix for the rename labels of the samples (e.g. Treated_1 if default)")
    parser.add_argument("-O", "--output", type=int, default=f"{os.getcwd()}", help="Directory to save the generated graphs")
    parser.add_argument("-v", "--save-mtx", action="store_true", help="Saves the renamed and normalized count matrix")
    parser.add_argument("-i", "--index-col", type=int, default=0, help="Column that will become the row index")
    parser.add_argument("-p", "--palette-colors", type=str, default="Paired", help="Seaborn palette colors")
    parser.add_argument("--MDS", action="store_true", help="Plots a metric-MDS among samples")
    parser.add_argument("--PLS", action="store_true", help="Plots a PLS Regression of the gene expression data")

    return parser.parse_args()

def get_count_mtx(count_path: str, sep: int = '\t', index_col: int = 0) -> DataFrame:
    if not os.path.isfile(count_path):
        lg.error(f"[GET] File {os.path.basename(count_path)} not available at {os.path.dirname(count_path)}")
        sys.exit(1)

    count = read_csv(count_path, sep=sep, index_col=index_col)
    if count.index.name:
        count.index.name = None
    
    return count

# ============================================================================= #
#                           STANDARIZING COUNT MATRIX
# ============================================================================= #

def __standarize_df(
        count_path: str, 
        sep: int = '\t', 
        index_col: int = 0,
        order: str = '4-6,1-3,7',
        conditions: str = 'Treated,Untreated',
        blocks: str = '3,4',
        start_at: int = 1,
) -> DataFrame:
    count = get_count_mtx(count_path=count_path, sep=sep, index_col=index_col)
    count.index.name = None
    lg.info("[STANDARIZE] Count matrix succesfully charged in a DataFrame")

    try:

        idx = []
        for e in order.split(',').strip():
            if re.fullmatch(r"(\d+)-(\d+)", e):
                a, b = e.split('-').strip()
                idx.extend([x for x in range(int(a), int(b)+1)])
            else:
                idx.append(int(e))

    except Exception as e:
        lg.error("[STANDARIZE] No valid indexes to reorder the DataFrame")

    try:
        count = count.iloc[:, idx]
        labels = []
        blocks = blocks.split(",").strip()
        for i, c in enumerate(conditions.split(",").strip()):
            labels.extend([f"{c}_{i}" for i in range(start_at, blocks[i]+start_at)])

        count.columns = labels

    except Exception as e:
        lg.error("[STANDARIZE] Fail to reorganize the DataFrames with the given index")

    return count

def __normalize_mtx(count_stand: DataFrame, logbase: int = 2) -> DataFrame:
    try:
        count_log2 = np.log2(count_stand + 1.0)
    except Exception as e:
        lg.error(f"[NORMALIZE] Fail to apply logarithmic (base {logbase}) normalization")

    return count_log2

def __melt_mtx(count_log2: DataFrame) -> DataFrame:
    try:
        count_log2_melt = count_log2.melt(var_name="Sample", value_name="Value")
        count_log2_melt["Condition"] = count_log2_melt['Sample'].str.replace('_.$','', regex=True)
    except Exception as e:
        lg.error("[MELT] Fail to melt samples and log values")

    lg.info("[MELT] The DataFrame was succesfully melted")
    return count_log2_melt

# ============================================================================= #
#                                    CALC
# ============================================================================= #

def fit_pca(expression: DataFrame) -> tuple[PCA, np.ndarray]:
    pca = PCA(n_components=2, svd_solver="full", tol=1e-9, random_state=5)

    return pca, pca.fit_transform(expression.values) # Return model and coords

def corrp(m: np.ndarray, axis: int = 1) -> np.ndarray:
    Z = zscore(m, axis=axis, ddof=1)
    if np.isnan(Z):
        msg = "[CORRP] NaNs found after Zscore (constant samples?)"
        lg.error(msg=msg)
        raise ValueError(msg)
    
    p = m.shape[1]
    R = (Z @ Z.T) / (p-1) # Pearson correlation based on thee covariance
    np.fill_diagonal(R, 1.0)
    R = np.clip(R, -1.0, 1.0)
    D = 1 - R # Distance based on Pearson correlation
    return (D + D.T) / 2 # Force symmetry

def fit_mds(expression_log2: DataFrame, metric: str = "correlation", axis: int = 1, workers: int = 4) -> tuple[MDS, np.ndarray]:
    valid_metrics = {"correlation", "manhattan", "euclidean"}

    if not (metric := metric.lower()) in valid_metrics:
        valid = "".join(valid_metrics)
        msg = f"[FIT_MDS] No valid metric was provided, must be {valid}"
        lg.error(msg=msg)
        raise ValueError(msg)
    
    m = expression_log2.values
    conditions, genes = m.shape

    if genes > conditions and metric == "correlation":
        D = corrp(m=m, axis=1)
    else:
        total = mp.cpu_count()
        v = pdist(m, metric=metric, workers=max(2, total // 16))
        D = squareform(v)

    mds = MDS(
        n_components=2, 
        dissimilarity='precomputed', 
        metric=True,
        n_init=15, 
        max_iter=1000, 
        eps=1e-08, 
        normalized_stress='auto',
        random_state=5
    )
    
    return mds, mds.fit_transform(D) # Return model and coords

# ============================================================================= #
#                                     PLOT
# ============================================================================= #

def boxplot(
        data: DataFrame,
        x: str = "Sample",
        y: str = "Value",
        hue: str = "Condition",
        title: str = "Distribution per condition",
        color_palette: str = "Paired",
        rotation: bool = True,
        svfig: bool = True,
        dpi: int = 300,
        filename: str = "condition_distributions.png",
        format: str = "png,svg",
        outdir: str = "results/Images"
        ) -> None:
    valid_formats = {"png", "pdf", "svg", "eps", "ps"}
    palette = sns.color_palette(color_palette)
    sns.boxplot(data=data, x=x, y=y, hue=hue, palette=[palette[i] for i, _ in enumerate(set(data[hue]))])
    plt.tight_layout()
    plt.title(title)
    plt.legend(
        bbox_to_anchor=(1.02, 0.5),   # Position outside the axe
        loc='upper left',           # Legend location/position
        borderaxespad=0,
        fontsize = 8,
        title_fontsize=10,
        labelcolor = "black"
    )
    if rotation:
        plt.xticks(rotation=0.45)
    if svfig:
        try:
            formats = {}
            if (fmt := filename.split(".")[-1]) and fmt in valid_formats:
                formats.update(fmt)
            if (fmt := format.split(",").strip()):
                formats.update(f for f in fmt if f in valid_formats)

        except Exception as e:
            lg.error("[BOXPLOT] The saving format couldn't be applied")
            
        for fmt in formats:
            plt.savefig(f'{outdir}/{filename}.{fmt}', dpi=dpi, bbox_inches='tight')

        lg.info(f"[BOXPLOT] The plots were saved in {outdir} with the following formats: {formats}")
    else:
        plt.show()
        lg.info("[BOXPLOT] The plots were displayed in the user stout but not saved")

    return

def data_kde(
        data: DataFrame,
        x: str = "Value",
        hue: str = "Sample",
        col: str = "Condition",
        fill: bool = True,
        alpha: float = 0.25,
        suptitle: str = "Distribution per condition",
        color_palette: str = "Paired",
        svfig: bool = True,
        dpi: int = 300,
        filename: str = "condition_distributions.png",
        format: str = "png,svg",
        outdir: str = "results/Images"
) -> None:
    valid_formats = {"png", "pdf", "svg", "eps", "ps"}
    palette = sns.color_palette(color_palette)
    sns.displot(data=data, x=x, hue=hue, col=col, fill=fill, alpha=alpha, palette=color_palette)
    plt.suptitle(suptitle)
    plt.tight_layout()

    if svfig:
        try:
            formats = {}
            if (fmt := filename.split(".")[-1]) and fmt in valid_formats:
                formats.update(fmt)
            if (fmt := format.split(",").strip()):
                formats.update(f for f in fmt if f in valid_formats)

        except Exception as e:
            lg.error("[KDE] The saving format couldn't be applied")
            
        for fmt in formats:
            plt.savefig(f'{outdir}/{filename}.{fmt}', dpi=dpi, bbox_inches='tight')

        lg.info(f"[KDE] The plots were saved in {outdir} with the following formats: {formats}")
    else:
        plt.show()
        lg.info("[KDE] The plots were displayed in the user stout but not saved")

    return

# ============================================================================= #
#                                     MAIN
# ============================================================================= #

def main():
    args = my_parser()

if __name__ == "__main__":
    main()