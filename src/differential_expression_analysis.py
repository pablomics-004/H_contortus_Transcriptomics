"""
Differential Expression Analysis
--------------------------------
Receives a raw count matrix corrected by pyCombat and uses the DESeq2 algorithm to capture the
differential expression analysis.

Author: Pablo Salazar-Méndez
Version: 1.0.0
Date: 2025-11-23

This script seeks to automate the Differential Expression Analysis of raw count matrix by having
the next returns:

    - Genes classficiation and their associated statistics matrix
    - Differentially expressed genes IDs (UP, DOWN) table
    - Volcano plots
    - Heatmap

Dependencies:
    - matplotlib.pyplot
    - seaborn
    - pandas
    - numpy
    - pydeseq2
    - visualization
    - scipy or fastcluster

The latter dependency is located within this same directory in the github page of this pipeline.

Example usage:
    python3 src/differential_expression_analysis.py \
        -f results/filtered_counts.tsv \
        -D results/metadata.tsv -O results/Images/ \
        --DE_matrix_name diffexp_matrix_no_annot_lfc1.tsv \
        -L 1 -P 1e-2 -t 3 -S \
        --bar_name barplot_diffexp_no_annot_lfc1 \
        --volcano_name vulcano_diffexp_ids_lfc1 \
        --image_formats jpg --dpi 300 -p Spectral
"""

from argparse import Namespace, ArgumentParser
from visualization import get_count_mtx
from contextlib import redirect_stdout
from pydeseq2.dds import DeseqDataSet
from pandas import DataFrame, Series
from pydeseq2.ds import DeseqStats
from datetime import datetime

import matplotlib.pyplot as plt
import seaborn as sns
import logging as lg
import numpy as np
import os
import io

# ============================================================================= #
#                                   UTILITY
# ============================================================================= #

logger = lg.getLogger("differential_expression")

def setup_logger(level=lg.INFO):
    if logger.handlers:
        return logger

    logger.setLevel(level)
    logger.propagate = False

    log_filename = f"differential_expression_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    fh = lg.FileHandler(log_filename)
    sh = lg.StreamHandler()

    formatter = lg.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    sh.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(sh)

    return logger

def my_parser() -> Namespace:
    parser = ArgumentParser(description="Differential Expression Analysis pipeline")
    parser.add_argument("-f", "--filtered_counts", type=str, required=True, help="Filtered raw count standarize matrix for DEA")
    parser.add_argument("-D", "--design_matrix", type=str, help="Design matrix location within directories")
    parser.add_argument("--build_design_matrix", action="store_true", help="If the design matrix is not available, the pipeline builds it")
    parser.add_argument("-L", "--log2foldchange", type=int, default=2, help="LFC criteria to determine up- and downregulated genes")
    parser.add_argument("-P", "--padjust", type=float, default=1e-2, help="padjust value criteria for statistical significance")
    parser.add_argument("-O", "--outdir", type=str, default=f"{os.getcwd()}", help="Stores the directory for the results")
    parser.add_argument("--DE_matrix_name", type=str, default="diffexp_matrix.tsv", help="Name for the differential expression analysis matrix")
    parser.add_argument("-s", "--separator", type=str, default="\t", help="Separator of the filtered counts matrix. The output file of the pipeline will have the same separator")
    parser.add_argument("-i", "--index_col", type=int, default=0, help="Column that will become the row index")
    parser.add_argument("-S", "--save_ids", action="store_true", help="Saves a DataFrame with the gene_ids of the DE genes")
    parser.add_argument("--bar_name", type=str, default="bar_genes_classification.jpg", help="Bar plot name")
    parser.add_argument("--volcano_name", type=str, default="vulcano_diffexp.jpg", help="Vulcano plot name")
    parser.add_argument("-t", "--top_genes", type=int, default=0, help="Top differentialy expressed genes")
    parser.add_argument("--gene_name_col", type=str, default="", help="Name of the gene name column associated with the gene ID")
    parser.add_argument("--heatmap_name", type=str, default="heatmap_diffexp.jpg", help="Heatmap plot name")
    parser.add_argument("--image_formats", type=str, default="jpg,svg", help="Format for saving the plots separated by a coma ','")
    parser.add_argument("--dpi", type=int, default=300, help="Resolution of the image in dots per inch (dpi)")
    parser.add_argument("-p", "--palette_colors", type=str, default="Spectral", help="Seaborn palette colors")


    return parser.parse_args()

def build_design_mtx(stand_count_matrix: DataFrame) -> DataFrame:
    idx = stand_count_matrix.columns
    conditions = idx.str.replace("_[0-9]+$", "", regex=True)
    return DataFrame(conditions, columns=["condition"], index=idx)

def get_label(df: DataFrame, name_col: str = "gene_name") -> Series:
    return df[name_col].astype(str) if name_col and name_col in df.columns else df.index.astype(str)

# ============================================================================= #
#                                   pyDESEQ2
# ============================================================================= #

def diffexp(
        counts: DataFrame, 
        design_matrix: DataFrame, 
        logfold: int = 1, 
        padj: float = 1e-2
) -> tuple[DataFrame, DataFrame]:
    counts_transpose = counts.T
    conditions = list(design_matrix.iloc[:,0].unique())
    reference = conditions[0]
    interest = conditions[1]
    lfc = "log2FoldChange"
    var = "condition"
    p = "padj"

    logger.info(f"[DIFFEXP] Design matrix for pyDESeq2:\n{design_matrix}")

    with redirect_stdout(io.StringIO()) as f:
        dds = DeseqDataSet(
            counts=counts_transpose,
            metadata=design_matrix,
            design=f"~ {var}"
        )
        dds.deseq2()
        stats = DeseqStats(dds=dds, contrast=[var, interest, reference])
        stats.summary()

    logger.info(f"[DIFFEXP] Stats summary captured:\n{f.getvalue()}")

    results = stats.results_df

    try:
        cond = [
            (results[p] <= padj) & (results[lfc] >= logfold) & (results[p].notna()),
            (results[p] <= padj) & (results[lfc] <= -logfold) & (results[p].notna())
        ]
        opc = ("UP", "DOWN")
        results["Expression"] = np.select(cond, opc, default="NON-DE")
        results.loc[results[p].isna(), "Expression"] = "NO-PADJ"
    except Exception as e:
        msg = f"[DIFFEXP] Fail to classify expression due to the following error: {e}"
        logger.error(msg=msg)
        raise ValueError(msg)
    
    logger.info("[DIFFEXP] Succesfully applied DESeq2 algorithm")
    
    return results, dds.layers['normed_counts']

def annotate_top_genes(ax, stats: DataFrame, n: int = 20, name_col: str = "") -> None:
    up = stats[stats["Expression"] == "UP"].copy()
    down = stats[stats["Expression"] == "DOWN"].copy()

    up = up.sort_values(["padj", "log2FoldChange"], ascending=[True, False]).head(n)
    down = down.sort_values(["padj", "log2FoldChange"], ascending=[True, True]).head(n)

    up_labels = get_label(df=up, name_col=name_col)
    down_labels = get_label(df=down, name_col=name_col)
    batches = [(up, up_labels, 0.12, 0.06, "left"),(down, down_labels, -0.12, 0.06, "right")]

    for df, df_labels, dx, dy, ha in batches:
        for (gene, row), label in zip(df.iterrows(), df_labels):
            x, y = row["log2FoldChange"], row["log10Neg"]
            ax.annotate(
                label,
                xy=(x, y),
                xytext=(x + dx, y + dy),
                textcoords='data',
                arrowprops=dict(arrowstyle='-',lw=0.4,alpha=0.5),
                fontsize=7,
                ha=ha,
                va="center"
            )

# ============================================================================= #
#                                     PLOTS
# ============================================================================= #

def bar_freq(
        frequencies: DataFrame,
        svfig: bool = False,
        filename: str = "genes_classification",
        outdir: str = "results/",
        outfmt: str = "jpg,svg",
        color_palette: str = "Paired",
        dpi: int = 300
) -> None:
    valid_formats = {"png", "pdf", "svg", "eps", "ps", "jpg", "jpeg"}
    fig, ax = plt.subplots(figsize=(4, 6))
    sns.barplot(x=frequencies.index, y=frequencies.values, palette=color_palette, ax=ax)

    ax.set_yscale("log")
    ax.set_ylabel("Count (log scale)")
    ax.set_xlabel("Expression category")
    ax.set_title("Differential Expression summary")

    for i, value in enumerate(frequencies.values):
        ax.text(i, value * 1.05, str(value), ha='center', va='bottom', fontsize=11)

    fig.tight_layout()

    if svfig:
        try:
            formats = set()
            if (fmt := filename.split(".")[-1].strip()) and fmt in valid_formats:
                formats.add(fmt)
            if (fmt := outfmt.split(",")):
                formats.update(e for f in fmt if (e:=f.strip()) in valid_formats)

        except Exception as e:
            msg = "[BAR_FREQ] The saving format couldn't be applied"
            logger.error(msg)
            raise ValueError(msg)
            
        for fmt in formats:
            fname, _ = os.path.splitext(os.path.basename(filename))
            plt.savefig(f'{outdir}/{fname}.{fmt}', dpi=dpi, bbox_inches='tight')

        logger.info(f"[BAR_FREQ] The plot was saved in {outdir} with the following formats: {formats}")
    else:
        plt.show()
        logger.info("[BAR_FREQ] The plot was displayed in the user stout but not saved")

    plt.close(fig=fig)

    return

def volcano_plot(
    stats: DataFrame,
    lfc_criteria: int = 1,
    padj_criteria: float = 1e-2,
    gene_col_name: str = "gene_name",
    top_genes: int = 0,
    svfig: bool = False,
    filename: str = "genes_classification",
    outdir: str = "results/",
    outfmt: str = "jpg,svg",
    color_palette: str = "Spectral",
    dpi: int = 300
) -> None:
    valid_formats = {"png", "pdf", "svg", "eps", "ps", "jpg", "jpeg"}
    cmap = sns.color_palette(color_palette, as_cmap=True)
    log10_criteria = -np.log10(padj_criteria)
    order = ["NON-DE", "DOWN", "UP"]
    colors = {"UP" : cmap(1.0), "NON-DE" : (0.75, 0.75, 0.75, 0.3), "DOWN" : cmap(0.0)}
    fig, ax = plt.subplots(figsize=(10, 6))

    sns.scatterplot(
        data=stats,
        x="log2FoldChange",
        y="log10Neg",
        hue="Expression",
        hue_order=order,
        palette=colors,
        s=20,
        ax=ax
    )

    ax.set_xlabel("log2(Fold Change)")
    ax.set_ylabel("-log10(padj)")
    ax.set_title("Volcano Plot")
    ax.grid(True, alpha=0.5, linestyle="--")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axvline(lfc_criteria,  color='gray', linestyle='--', alpha=0.7)
    ax.axvline(-lfc_criteria, color='gray', linestyle='--', alpha=0.7)
    ax.axhline(log10_criteria, color='gray', linestyle='--', alpha=0.7)

    if top_genes:
        annotate_top_genes(ax=ax,stats=stats,n=top_genes,name_col=gene_col_name)

    fig.tight_layout()

    if svfig:
        try:
            formats = set()
            if (fmt := filename.split(".")[-1].strip()) and fmt in valid_formats:
                formats.add(fmt)
            if (fmt := outfmt.split(",")):
                formats.update(e for f in fmt if (e:=f.strip()) in valid_formats)

        except Exception as e:
            msg = "[VOLCANO] The saving format couldn't be applied"
            logger.error(msg)
            raise ValueError(msg)
            
        for fmt in formats:
            fname, _ = os.path.splitext(os.path.basename(filename))
            plt.savefig(f'{outdir}/{fname}.{fmt}', dpi=dpi, bbox_inches='tight')

        logger.info(f"[VOLCANO] The plot was saved in {outdir} with the following formats: {formats}")
    else:
        plt.show()
        logger.info("[VOLCANO] The plot was displayed in the user stout but not saved")

    plt.close(fig=fig)

    return

def heatmap(
        data: DataFrame,
        zscore: int = 0,
        svfig: bool = False,
        filename: str = "genes_classification",
        outdir: str = "results/",
        outfmt: str = "jpg,svg",
        color_palette: str = "Spectral",
        dpi: int = 300
) -> None:
    valid_formats = {"png", "pdf", "svg", "eps", "ps", "jpg", "jpeg"}
    cmap = sns.color_palette(color_palette, as_cmap=True)

    cg = sns.clustermap(
        data,
        cmap=cmap,
        z_score=zscore,
        row_cluster=True,
        col_cluster=True,
        figsize=(10, 8),
        center=0
    )

    cg.figure.suptitle("Heatmap of Differentially Expressed Genes", y=1.02)
    cg.cax.set_ylabel(
        "Z-score (log1p counts)\n(Down  ←   →  Up)",
        rotation=90,
        ha="center",
        va="center",
        labelpad=12
    )
    cg.figure.tight_layout()
    cg.figure.subplots_adjust(top=0.95)

    if svfig:
        try:
            formats = set()
            if (fmt := filename.split(".")[-1].strip()) and fmt in valid_formats:
                formats.add(fmt)
            if (fmt := outfmt.split(",")):
                formats.update(e for f in fmt if (e:=f.strip()) in valid_formats)

        except Exception as e:
            msg = "[HEATMAP] The saving format couldn't be applied"
            logger.error(msg)
            raise ValueError(msg)
            
        for fmt in formats:
            fname, _ = os.path.splitext(os.path.basename(filename))
            cg.figure.savefig(f"{outdir}/{fname}.{fmt}", dpi=dpi, bbox_inches="tight")

        logger.info(f"[HEATMAP] The plot was saved in {outdir} with the following formats: {formats}")
    else:
        plt.show()
        logger.info("[HEATMAP] The plot was displayed in the user stout but not saved")

    plt.close(cg.figure)

    return

# ============================================================================= #
#                                     MAIN
# ============================================================================= #

def main():
    setup_logger()
    args = my_parser()

    os.makedirs(args.outdir, exist_ok=True)

    logger.info("[MAIN] Starting pipeline")
    counts = get_count_mtx(args.filtered_counts, sep=args.separator, index_col=args.index_col)
    counts = counts.round(0).clip(lower=0).astype(int)
    logger.info(f"[MAIN] Filtered count matrix fetched from {args.filtered_counts}")

    rows, cols = counts.shape
    if cols > rows:
        msg = "Your count matrix has more columns than rows  — are samples columns?"
        logger.error(msg=msg)
        raise ValueError(msg)

    if args.build_design_matrix:
        metadata = build_design_mtx(counts)
    else:
        metadata = get_count_mtx(args.design_matrix, sep=args.separator, index_col=args.index_col)
    metadata.columns = ["condition"]
    logger.info(f"[MAIN] Successfully obtained design matrix")

    # ============================================ pyDESeq2 ============================================ #

    logger.info("[MAIN] Initiating differentia expression analysis with pyDESeq2")
    dea, norm_counts = diffexp(counts=counts, design_matrix=metadata, logfold=args.log2foldchange, padj=args.padjust)
    dea["log10Neg"] = -np.log10(dea["padj"].clip(lower=1e-300))
    freqs = dea["Expression"].value_counts()
    logger.info("[MAIN] Differential expression processing finished")

    dea_clean = dea.dropna(subset=["padj"])
    norm_counts_df = DataFrame( # samples x genes
        norm_counts,
        index=counts.columns, 
        columns=counts.index
    )

    if args.separator in {"\t", r"\t"}:
        sep, ext = "\t", "tsv"
    elif args.separator == ",":
        sep, ext = ",", "csv"
    else:
        sep, ext = "\t", "tsv"
        logger.error("[MAIN] Unsupported separator, saving as TSV")
    
    parent = os.path.dirname(os.path.normpath(args.outdir)) or os.path.normpath(args.outdir)
    fname = os.path.splitext(os.path.basename(args.DE_matrix_name))[0]
    dea_clean.to_csv(os.path.join(parent, f"{fname}.{ext}"), sep=sep)
    logger.info(f"[MAIN] Differential expression stats table saved at {parent}")

    updown_genes = dea_clean[dea_clean["Expression"].isin(["UP","DOWN"])]

    if args.save_ids:
        diff_exp_genes = (
            updown_genes
            .reset_index()
            .rename(columns={"index" : "gene_id"})
            [["gene_id", "Expression"]]
        )
        diff_exp_genes.to_csv(os.path.join(parent, f"{fname}_only_diffexp.{ext}"), sep=sep)
        logger.info(f"[MAIN] Matrix for differential expression IDs saved at {parent}")

    # ============================================== PLOTS ============================================== #
    
    bar_freq(
        frequencies=freqs, 
        svfig=True, 
        filename=args.bar_name, 
        outdir=args.outdir,
        outfmt=args.image_formats,
        color_palette=args.palette_colors,
        dpi=args.dpi
    )

    # Volcano
    volcano_plot(
        stats=dea_clean,
        lfc_criteria=args.log2foldchange,
        padj_criteria=args.padjust,
        gene_col_name=args.gene_name_col,
        top_genes=args.top_genes,
        svfig=True,
        filename=args.volcano_name,
        outdir=args.outdir,
        outfmt=args.image_formats,
        color_palette=args.palette_colors,
        dpi=args.dpi
    )

    # Heatmap
    gene_names = updown_genes.index
    heatmap_data = np.log1p(norm_counts_df.loc[:, gene_names]).T

    heatmap(
        data=heatmap_data,
        zscore=0,
        svfig=True,
        filename=args.heatmap_name,
        outdir=args.outdir,
        outfmt=args.image_formats,
        color_palette=args.palette_colors,
        dpi=args.dpi
    )

    logger.info("[MAIN] Differential Expression Analysis pipeline finished")
    return

if __name__ == "__main__":
    main()