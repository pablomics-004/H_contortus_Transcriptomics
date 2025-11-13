#!/usr/bin/env python3
'''
Data visualization before analysis
----------------------------------
Given a count matrix from gene expresion data, allows the visualization of the data before
doing a Differential Gene Expression analysis. It supose the next structure of the count matrix:

gene_id,gene_name,sample1,sample2,...,samplen

Authors: Yael Montiel-Vargas & Pablo Salazar-Méndez
Version: 1.0.0
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
    parser.add_argument("-s", "--separator", type=str, default="\t", help="Separator of the count matrix (default is a tab)")
    parser.add_argument("-l", "--log_base", type=int, default=2, help="Base for the logarithmic normalization")
    parser.add_argument("-o", "--order_col", required=True, type=str, help="Final order of the columns within the DataFrame (Pandas)")
    parser.add_argument("-f", "--format_images", type=str, default="png,svg", help="Format (or formats) for saving the images")
    parser.add_argument("--start_at", type=int, default=1, help="Sufix for the rename labels of the samples (e.g. Treated_1 if default)")
    parser.add_argument("-O", "--output", type=str, default=f"{os.getcwd()}", help="Directory to save the generated graphs")
    parser.add_argument("-d","--dpi", type=int, default=300, help="Resolution of the image in dots per inch (dpi)")
    parser.add_argument("-v", "--save_mtx", action="store_true", help="Saves the renamed and normalized count matrix")
    parser.add_argument("-i", "--index_col", type=int, default=0, help="Column that will become the row index")
    parser.add_argument("-p", "--palette_colors", type=str, default="Paired", help="Seaborn palette colors")
    parser.add_argument("--mds", action="store_true", help="Plots a metric-MDS among samples")
    parser.add_argument("--metric", type=str, default="correlation", help="Metric to compute the distance matrix for the MDS among samples")
    parser.add_argument("-w", "--workers", type=int, default=4, help="Number of worker to compute the distance matrix")
    parser.add_argument("--boxplot_name", type=str, default="distribution_per_condition.png", help="Name for box plot image")
    parser.add_argument("--kde_name", type=str, default="kde_per_condition.png", help="Name for KDE plot image")
    parser.add_argument("--pca_name", type=str, default="samples_pca.png", help="Name for saving the samples PCA plot")
    parser.add_argument("--mds_name", type=str, default="samples_mds.png", help="Name for saving the samples MDS plot")

    return parser.parse_args()

def get_count_mtx(count_path: str, sep: str = "\t", index_col: int = 0) -> DataFrame:
    if not os.path.isfile(count_path):
        msg = f"[GET] File {os.path.basename(count_path)} not available at {os.path.dirname(count_path)}"
        lg.error(msg=msg)
        raise ValueError(msg)

    count = read_csv(count_path, sep=sep, index_col=index_col)
    if count.index.name:
        count.index.name = None
    
    return count

# ============================================================================= #
#                           STANDARIZING COUNT MATRIX
# ============================================================================= #

def __standarize_df(
        count_path: str, 
        sep: str = "\t", 
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
        for e in order.split(','):
            e = e.strip()
            if re.fullmatch(r"(\d+)-(\d+)", e):
                a, b = e.split('-')
                a = a.strip(); b = b.strip()
                idx.extend([x for x in range(int(a)-1, int(b))])
            else:
                idx.append(int(e)-1)

    except Exception as e:
        msg = "[STANDARIZE] No valid indexes to reorder the DataFrame"
        lg.error(msg=msg)
        raise KeyError(msg)

    try:
        count = count.iloc[:, idx]
        labels = []
        blocks = [int(b.strip()) for b in blocks.split(",")]
        for i, c in enumerate(conditions.split(",")):
            labels.extend([f"{c.strip()}_{j}" for j in range(start_at, blocks[i]+start_at)])

        count.columns = labels

    except Exception as e:
        msg = "[STANDARIZE] Fail to reorganize the DataFrames with the given index"
        lg.error(msg=msg)
        raise ValueError(msg)
    
    if sum(blocks) != len(labels):
        msg = "[STANDATIZE] Discrepancy of the number of blocks and reorganized columns"
        lg.error(msg=msg)
        raise ValueError(msg)

    return count

def __normalize_mtx(count_stand: DataFrame, logbase: int = 2) -> DataFrame:
    if logbase < 0:
        msg = f"[NORMALIZE] Negative logarithm bases aren't supported"
        lg.error(msg=msg)
        raise ValueError(msg)
    try:      

        match logbase:
            case 0:
                count_log = np.log(count_stand + 1.0)
            case 2:
                count_log = np.log2(count_stand + 1.0)
            case 10:
                count_log = np.log10(count_stand + 1.0)
            case _:
                count_log = np.log(count_stand + 1.0) / np.log(logbase)

    except Exception as e:
        msg = f"[NORMALIZE] Fail to apply logarithmic (base {logbase}) normalization"
        lg.error(msg=msg)
        raise ValueError(msg)

    return count_log

def __melt_mtx(count_log2: DataFrame) -> DataFrame:
    try:
        count_log2_melt = count_log2.melt(var_name="Sample", value_name="Value")
        count_log2_melt["Condition"] = count_log2_melt['Sample'].str.replace('_.$','', regex=True)
    except Exception as e:
        msg = "[MELT] Fail to melt samples and log values"
        lg.error(msg=msg)
        raise ValueError(msg)

    lg.info("[MELT] The DataFrame was succesfully melted")
    return count_log2_melt

# ============================================================================= #
#                                    CALC
# ============================================================================= #

def fit_pca(expression: DataFrame) -> tuple[np.ndarray, np.ndarray]:
    pca = PCA(n_components=2, svd_solver="full", tol=1e-9, random_state=5)
    
    coords = pca.fit_transform(expression.values)
    exp_var = pca.explained_variance_ratio_ * 100

    return exp_var, coords # Return explained variance and coords

def corrp(m: np.ndarray, axis: int = 1) -> np.ndarray:
    Z = zscore(m, axis=axis, ddof=1)
    if np.isnan(Z).any():
        msg = "[CORRP] NaNs found after zscore (constant samples?)"
        lg.error(msg=msg)
        raise ValueError(msg)
    
    p = m.shape[1]
    R = (Z @ Z.T) / (p-1) # Pearson correlation based on the sample covariance
    np.fill_diagonal(R, 1.0)
    R = np.clip(R, -1.0, 1.0)
    D = 1 - R # Distance based on Pearson correlation
    return (D + D.T) / 2 # Force symmetry

def fit_mds(expression_log2: DataFrame, metric: str = "correlation", workers: int = 4) -> tuple[float, np.ndarray]:
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
        v = pdist(m, metric=metric, workers=workers)
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

    coords = mds.fit_transform(D)
    coords = coords - coords.mean(axis=0, keepdims=True)
    
    return mds.stress_, coords # Return stress value, coords and distance matrix

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
            formats = set()
            if (fmt := filename.split(".")[-1].strip()) and fmt in valid_formats:
                formats.add(fmt)
            if (fmt := format.split(",")):
                formats.update(e for f in fmt if (e:=f.strip()) in valid_formats)

        except Exception as e:
            msg = "[BOXPLOT] The saving format couldn't be applied"
            lg.error(msg=msg)
            raise ValueError(msg)
            
        for fmt in formats:
            fname, _ = os.path.splitext(os.path.basename(filename))
            plt.savefig(f'{outdir}/{fname}.{fmt}', dpi=dpi, bbox_inches='tight')

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
    sns.displot(data=data, x=x, hue=hue, col=col, fill=fill, alpha=alpha, palette=color_palette)
    plt.suptitle(suptitle)
    plt.tight_layout()

    if svfig:
        try:
            formats = set()
            if (fmt := filename.split(".")[-1].strip()) and fmt in valid_formats:
                formats.add(fmt)
            if (fmt := format.split(",")):
                formats.update(e for f in fmt if (e:=f.strip()) in valid_formats)

        except Exception as e:
            msg = "[KDE] The saving format couldn't be applied"
            lg.error(msg=msg)
            raise ValueError(msg)
            
        for fmt in formats:
            fname, _ = os.path.splitext(os.path.basename(filename))
            plt.savefig(f'{outdir}/{fname}.{fmt}', dpi=dpi, bbox_inches='tight')

        lg.info(f"[KDE] The plots were saved in {outdir} with the following formats: {formats}")
    else:
        plt.show()
        lg.info("[KDE] The plots were displayed in the user stout but not saved")

    return

def plot_mds(
        mds_df: DataFrame, 
        stress: float,
        x: str = "Dimension 1",
        y: str = "Dimension 2",
        hue: str = "Condition",
        style : str = "Condition",
        title: str = "MDS",
        edgecolor: str = "purple",
        color_palette: str = "Paired",
        svfig: bool = True,
        dpi: int = 300,
        filename: str = "samples_mds.png",
        format: str = "png,svg",
        outdir: str = "results/Images"
) -> None:
    valid_formats = {"png", "pdf", "svg", "eps", "ps"}
    sns.scatterplot(
        data=mds_df,
        x=x,
        y=y,
        hue="Condition",
        style="Condition",
        edgecolor="purple"
    )
    plt.legend(
        [f"stress = {stress}"],
        #bbox_to_anchor=(0.7, 0.7),   # Position outside the axe
        loc='upper right',           # Legend location/position
        borderaxespad=0,
        fontsize = 12,
        title_fontsize=10,
        labelcolor = "black"
    )
    plt.title(title)
    plt.tight_layout()
    plt.axis("equal")
    plt.grid(True, alpha=0.5, linestyle="--")

    if svfig:
        try:
            formats = set()
            if (fmt := filename.split(".")[-1].strip()) and fmt in valid_formats:
                formats.add(fmt)
            if (fmt := format.split(",")):
                formats.update(e for f in fmt if (e:=f.strip()) in valid_formats)

        except Exception as e:
            msg = "[MDS] The saving format couldn't be applied"
            lg.error(msg=msg)
            raise ValueError(msg)
            
        for fmt in formats:
            fname, _ = os.path.splitext(os.path.basename(filename))
            plt.savefig(f'{outdir}/{fname}.{fmt}', dpi=dpi, bbox_inches='tight')

        lg.info(f"[MDS] The plot was saved in {outdir} with the following formats: {formats}")
    else:
        plt.show()
        lg.info("[MDS] The plot was displayed in the user stout but not saved")

    return

def plot_pca(
        pca_df: DataFrame, 
        explained_variance: np.ndarray,
        x: str = "PC1",
        y: str = "PC2",
        hue: str = "Condition",
        style : str = "Condition",
        title: str = "PCA",
        edgecolor: str = "purple",
        color_palette: str = "Paired",
        svfig: bool = True,
        dpi: int = 300,
        filename: str = "samples_pca.png",
        format: str = "png,svg",
        outdir: str = "results/Images"
) -> None:
    valid_formats = {"png", "pdf", "svg", "eps", "ps"}
    sns.scatterplot(
        data=pca_df,
        x=f"{x} ({explained_variance[0]:.2f}%)",
        y=f"{y} ({explained_variance[1]:.2f}%)",
        hue="Condition",
        style="Condition",
        edgecolor="purple"
    )
    plt.legend(
        #bbox_to_anchor=(0.7, 0.7),   # Position outside the axe
        loc='upper right',           # Legend location/position
        borderaxespad=0,
        fontsize = 12,
        title_fontsize=10,
        labelcolor = "black"
    )
    plt.title(title)
    plt.tight_layout()
    plt.grid(True, alpha=0.5, linestyle="--")

    if svfig:
        try:
            formats = set()
            if (fmt := filename.split(".")[-1].strip()) and fmt in valid_formats:
                formats.add(fmt)
            if (fmt := format.split(",")):
                formats.update(e for f in fmt if (e:=f.strip()) in valid_formats)

        except Exception as e:
            msg = "[PCA] The saving format couldn't be applied"
            lg.error(msg)
            raise ValueError("[PCA] The saving format couldn't be applied")
            
        for fmt in formats:
            fname, _ = os.path.splitext(os.path.basename(filename))
            plt.savefig(f'{outdir}/{fname}.{fmt}', dpi=dpi, bbox_inches='tight')

        lg.info(f"[PCA] The plot was saved in {outdir} with the following formats: {formats}")
    else:
        plt.show()
        lg.info("[PCA] The plot was displayed in the user stout but not saved")

    return

# ============================================================================= #
#                                     MAIN
# ============================================================================= #

def main():
    args = my_parser()

    count_matrix = __standarize_df(args.matrix, args.separator, args.index_col, args.order_col, args.conditions, args.blocks, args.start_at)
    norm_count_mtx = __normalize_mtx(count_matrix, args.log_base)
    melt_mtx = __melt_mtx(norm_count_mtx)

    if args.save_mtx:
        path = os.path.dirname(args.matrix)
        fname, _ = os.path.splitext(os.path.basename(args.matrix))

        match args.separator:
            case "\t":
                fmt = ".tsv"
            case ",":
                fmt = ".csv"
            case _:
                msg = f"[MAIN] Invalid format for a table"
                lg.error(msg=msg)
                raise ValueError(msg)
            
        norm_count_mtx.to_csv(os.path.join(path, f"{fname}_log{args.log_base}_norm{fmt}"), sep=args.separator, index=False)
        lg.info(f"[SAVE MATRIX] Normalized (log{args.log_base}) count matrix saved at {path}")
    
    labels = melt_mtx.columns
    os.makedirs(args.output, exist_ok=True)

    boxplot(
        data=melt_mtx,
        x=labels[0], y=labels[1], hue=labels[2],
        color_palette=args.palette_colors,
        dpi=args.dpi,
        svfig=True if os.path.isdir(args.output) else False,
        filename=args.boxplot_name,
        format=args.format_images,
        outdir=args.output
    )

    data_kde(
        data=melt_mtx,
        x=labels[1], hue=labels[0], col=labels[2],
        suptitle=args.kde_name,
        color_palette=args.palette_colors,
        dpi=args.dpi,
        filename=args.kde_name,
        format=args.format_images,
        outdir=args.output
    )

    expression_transposed = norm_count_mtx.T
    conditions = norm_count_mtx.columns.str.replace('_.$', '', regex=True)

    # ============================================ PCA ============================================ #

    explain_var, principal_components = fit_pca(expression=expression_transposed)
    cols = [f'PC1 ({explain_var[0]:.2f}%)', f'PC2 ({explain_var[1]:.2f}%)']

    pca_df = DataFrame(data=principal_components, columns=cols, index=expression_transposed.index)
    pca_df["Condition"] = conditions

    plot_pca(
        pca_df=pca_df,
        explained_variance=explain_var,
        x="PC1", y="PC2",
        hue="Condition",
        style="Condition",
        title="Samples distribution",
        color_palette=args.palette_colors,
        dpi=args.dpi,
        filename=args.pca_name,
        format=args.format_images,
        outdir=args.output
    )

    if args.mds:
        sample_names = expression_transposed.index
        stress, mds_coords = fit_mds(expression_log2=expression_transposed, metric=args.metric, workers=args.workers)
        mds_df = DataFrame(mds_coords, index=sample_names, columns=sample_names, dtype=np.float64)

        plot_mds(
            mds_df=mds_df,
            stress=stress,
            x="dimension 1", y="dimension 2",
            hue="Condition",
            style="Condition",
            title="MDS",
            color_palette=args.palette_colors,
            dpi=args.dpi,
            filename=args.mds_name,
            format=args.format_images,
            outdir=args.output
        )


if __name__ == "__main__":
    main()