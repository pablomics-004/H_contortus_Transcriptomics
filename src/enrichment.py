#!/usr/bin/env python3
import os
import sys
import time
import logging
import argparse
import pandas as pd
import numpy as np
from Bio.KEGG import REST
from collections import Counter
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

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
    if not os.path.exists(path):
        logging.error(f"File not found: {path}")
        return False
    if os.path.getsize(path) == 0:
        logging.error(f"File is empty: {path}")
        return False
    return True

def validate_inputs(files):
    for f in files:
        if not file_exists_nonempty(f):
            logging.error(f"Aborting due to missing or empty file: {f}")
            sys.exit(1)
    logging.info("All input files validated.")

def load_inputs(args):
    logging.info("Loading input files...")
    ko_df = pd.read_csv(args.ko_file, sep="\t", names=['gene_id','KO'])
    background = set(pd.read_csv(args.background, header=None)[0])
    go_df = pd.read_csv(args.go_file, sep="\t", header=None, names=['gene_id','GO'])
    de_df = pd.read_csv(args.de_file, sep=",", header=0)
    return ko_df, background, go_df, de_df

def save_dataframe(df, outdir, filename, save_files, key):
    if key not in save_files:
        return

    # Crear outdir si no existe
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        path = os.path.join(outdir, filename)
    else:
        path = filename  # Guardar en directorio actual si no hay outdir

    df.to_csv(path, sep="\t", index=False)
    logging.info(f"Saved {key} to {path}")

# ============================
# GENE NAME MAPPING
# ============================
def get_gene_name_from_ko(ko_id: str) -> str:
    try:
        entry = REST.kegg_get(ko_id).read()
        for line in entry.split("\n"):
            if line.startswith("NAME"):
                return line.replace("NAME","").strip().split(";")[0]
        logging.warning(f"NAME not found for KO {ko_id}")
        return ko_id
    except Exception as e:
        logging.warning(f"Error retrieving KO {ko_id}: {e}")
        return ko_id

def map_gene_names_to_new_df(ko_df: pd.DataFrame, outdir, save_files):
    logging.info("Mapping gene names from KEGG...")
    gene_names = []
    for i, ko in enumerate(ko_df.get('KO', [])):
        gene_names.append(get_gene_name_from_ko(ko))
        time.sleep(0.2)
        if (i+1) % 50 == 0:
            logging.info(f"Processed {i+1}/{len(ko_df)} genes")
    names_df = pd.DataFrame({'gene_id': ko_df['gene_id'], 'KO': ko_df['KO'], 'gene_name': gene_names})
    save_dataframe(names_df, outdir, "gene_id_to_gene_name.tsv", save_files, "gene_names")
    return names_df

# ============================
# DE SPLIT
# ============================
def split_up_down(de_df, fc_thresholds=[1.0], alpha_values=[0.05], fc_col='log2FoldChange', pval_col='padj'):
    """
    Divide los genes DE en múltiples sets basados en rangos de FC y p-value.

    fc_thresholds: lista ordenada de thresholds, por ejemplo [1, 2]
    alpha_values: lista de p-values correspondientes a cada set, misma longitud que fc_thresholds
    """
    results = []
    de_sub = de_df.dropna(subset=[fc_col, pval_col])

    # Asegurarse de que thresholds están ordenados
    fc_thresholds_sorted = sorted(fc_thresholds)
    
    for i, (fc, alpha) in enumerate(zip(fc_thresholds_sorted, alpha_values)):
        if i == 0:
            # Primer set: genes entre 1er threshold y siguiente (si existe)
            upper_fc = fc_thresholds_sorted[i+1] if i+1 < len(fc_thresholds_sorted) else np.inf
            up = set(de_sub.loc[(de_sub[fc_col] >= fc) & (de_sub[fc_col] < upper_fc) & (de_sub[pval_col] < alpha), 'gene_id'])
            down = set(de_sub.loc[(de_sub[fc_col] <= -fc) & (de_sub[fc_col] > -upper_fc) & (de_sub[pval_col] < alpha), 'gene_id'])
        else:
            # Sets siguientes: genes mayores al threshold anterior
            lower_fc = fc_thresholds_sorted[i]
            up = set(de_sub.loc[(de_sub[fc_col] >= lower_fc) & (de_sub[pval_col] < alpha), 'gene_id'])
            down = set(de_sub.loc[(de_sub[fc_col] <= -lower_fc) & (de_sub[pval_col] < alpha), 'gene_id'])
        
        logging.info(f"Set {i+1}: FC >= {fc}, alpha={alpha} -> Up: {len(up)}, Down: {len(down)}")
        results.append({'fc': fc, 'alpha': alpha, 'up': up, 'down': down})
    
    return results


def extract_up_down_kos_gos(up_genes, down_genes, ko_df, go_df, outdir, save_files):
    up_kos = ko_df.loc[ko_df['gene_id'].isin(up_genes), 'KO'].dropna().unique()
    down_kos = ko_df.loc[ko_df['gene_id'].isin(down_genes), 'KO'].dropna().unique()
    up_gos = go_df.loc[go_df['gene_id'].isin(up_genes), 'GO'].dropna().unique()
    down_gos = go_df.loc[go_df['gene_id'].isin(down_genes), 'GO'].dropna().unique()

    if "KEGG_table" in save_files:
        save_dataframe(pd.Series(up_kos, name='KO'), outdir, "up_genes_KO.tsv", save_files, "KEGG_table")
        save_dataframe(pd.Series(down_kos, name='KO'), outdir, "down_genes_KO.tsv", save_files, "KEGG_table")
    if "GO_table" in save_files:
        save_dataframe(pd.Series(up_gos, name='GO'), outdir, "up_genes_GO.tsv", save_files, "GO_table")
        save_dataframe(pd.Series(down_gos, name='GO'), outdir, "down_genes_GO.tsv", save_files, "GO_table")
    return {'up_KO': up_kos, 'down_KO': down_kos, 'up_GO': up_gos, 'down_GO': down_gos}


def go_enrichment(genes_set: set, go_df: pd.DataFrame, de_df: pd.DataFrame, background: set, obo_file: str, prefix: str):
    """
    Realiza GO enrichment y genera una tabla extendida con genes y métricas DE.
    
    genes_set : set
        Genes de interés (up o down)
    go_df : pd.DataFrame
        DataFrame con columnas 'gene_id' y 'GO'
    de_df : pd.DataFrame
        DataFrame de resultados DE con 'gene_id', 'log2FoldChange', 'padj' y opcional 'gene_name'
    background : set
        Lista de genes de fondo
    obo_file : str
        Ruta al archivo GO OBO
    prefix : str
        Prefijo para guardar archivos
    """
    
    # 1. Preparar mapeo GO
    gene2go = go_df.groupby('gene_id')['GO'].apply(list).to_dict()
    obodag = GODag(obo_file)
    goea = GOEnrichmentStudy(background, gene2go, obodag, propagate_counts=True, alpha=0.05, methods=['fdr_bh'])
    
    # 2. Enriquecimiento GO
    results_all = goea.run_study(genes_set)
    results = []
    extended_rows = []
    
    for r in results_all:
        if r.p_fdr_bh < 0.05:
            results.append({
                'GO': r.GO,
                'name': r.name,
                'namespace': r.namespace,
                'p_fdr': r.p_fdr_bh,
                'ratio_in_study': r.ratio_in_study,
                'ratio_in_pop': r.ratio_in_pop,
                'n_genes': len(r.study_items)
            })
            
            # --- TABLA EXTENDIDA ---
            genes_in_term = [g for g in r.study_items if g in genes_set]
            for g in genes_in_term:
                if g in de_df['gene_id'].values:
                    row = de_df.loc[de_df['gene_id']==g].iloc[0]
                    extended_rows.append({
                        'GO': r.GO,
                        'GO_name': r.name,
                        'namespace': r.namespace,
                        'Gene_ID': g,
                        'Gene_Name': row.get('gene_name', g),
                        'log2FC': row.get('log2FoldChange', np.nan),
                        'pvalue': row.get('padj', np.nan),
                        'ratio_in_study': r.ratio_in_study,
                        'ratio_in_pop': r.ratio_in_pop
                    })
    
    results_df = pd.DataFrame(results)
    extended_df = pd.DataFrame(extended_rows)
    
    # 3. Guardar resultados
    results_df.to_csv(f"{prefix}_GO_enrichment.tsv", sep="\t", index=False)
    extended_df.to_csv(f"{prefix}_GO_extended_table.tsv", sep="\t", index=False)
    
    logging.info(f"{len(results_df)} significant GO terms saved for {prefix}")
    logging.info(f"Extended GO table with genes saved: {len(extended_df)} rows")
    
    # 4. Graficar
    if len(results_df) > 0:
        plot_go_enrichment(results_df, prefix)
    else:
        logging.warning(f"No significant GO terms found for {prefix}")
    
    return results_df, extended_df


def plot_go_enrichment(results: pd.DataFrame, prefix: str, top_n: int = 20):
    """Función mejorada para visualización de resultados GO
    
    top_n : int
        Número de términos GO a mostrar en las gráficas top
    """
    if len(results) == 0:
        return
    
    # Ordenar por p_fdr y tomar los top_n
    top_results = results.sort_values('p_fdr').head(top_n).copy()
    top_results['neg_log10_fdr'] = -np.log10(top_results['p_fdr'])
    top_results['enrichment_ratio'] = top_results['ratio_in_study'] / top_results['ratio_in_pop']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # --- GRÁFICA 1: Significancia estadística ---
    colors = plt.cm.Set2(np.linspace(0, 1, len(top_results)))
    bars1 = ax1.barh(range(len(top_results)), top_results['neg_log10_fdr'], color=colors)
    ax1.set_yticks(range(len(top_results)))
    ax1.set_yticklabels([name[:60] + '...' if len(name) > 60 else name for name in top_results['name']])
    ax1.set_xlabel('-log10(FDR)')
    ax1.set_title(f'Top {top_n} GO Terms - Statistical Significance\n({prefix})')
    
    for i, bar in enumerate(bars1):
        width = bar.get_width()
        n_genes = top_results.iloc[i].get('n_genes', np.nan)  # opcional si agregaste conteo
        ax1.text(width + 0.05, bar.get_y() + bar.get_height()/2, 
                 f'FDR={top_results.iloc[i]["p_fdr"]:.2e}\nn={n_genes}', 
                 ha='left', va='center', fontsize=8)
    
    # --- GRÁFICA 2: Ratio de enriquecimiento ---
    scatter = ax2.scatter(top_results['enrichment_ratio'], range(len(top_results)), 
                          c=top_results['neg_log10_fdr'], s=100, cmap='viridis', alpha=0.8)
    ax2.set_yticks(range(len(top_results)))
    ax2.set_yticklabels([])
    ax2.set_xlabel('Enrichment Ratio (Study/Population)')
    ax2.set_title('Enrichment Ratio vs GO Terms')
    ax2.axvline(x=1, color='red', linestyle='--', alpha=0.5, label='No enrichment')
    plt.colorbar(scatter, ax=ax2, label='-log10(FDR)')
    
    # Color por namespace
    namespace_colors = {'biological_process': '#1f77b4', 
                        'molecular_function': '#ff7f0e', 
                        'cellular_component': '#2ca02c'}
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                  markerfacecolor=color, markersize=8, label=ns.replace("_", " ").title())
                       for ns, color in namespace_colors.items()]
    ax2.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f"{prefix}_GO_enrichment_improved.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # --- Gráfica por namespace ---
    plot_go_by_namespace(results, prefix, top_n=top_n)


def plot_go_by_namespace(results: pd.DataFrame, prefix: str, top_n: int = 10):
    """Gráfica separada por categorías de GO"""
    if len(results) == 0:
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    namespaces = ['biological_process', 'molecular_function', 'cellular_component']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, (ns, color) in enumerate(zip(namespaces, colors)):
        ns_data = results[results['namespace'] == ns].sort_values('p_fdr').head(top_n).copy()
        if len(ns_data) > 0:
            ns_data['neg_log10_fdr'] = -np.log10(ns_data['p_fdr'])
            
            axes[i].barh(range(len(ns_data)), ns_data['neg_log10_fdr'], color=color)
            axes[i].set_yticks(range(len(ns_data)))
            axes[i].set_yticklabels([name[:50] + '...' if len(name) > 50 else name 
                                     for name in ns_data['name']], fontsize=9)
            axes[i].set_xlabel('-log10(FDR)')
            axes[i].set_title(f'{ns.replace("_", " ").title()}\n(n={len(ns_data)})')
    
    plt.tight_layout()
    plt.savefig(f"{prefix}_GO_by_namespace.png", dpi=300, bbox_inches='tight')
    plt.close()

def get_ko_pathway_mapping(ko_list):
    """
    Map KO IDs to pathways via KEGG REST API.
    Returns:
        ko2path: dict {ko_id: [pathway_id, ...]}
        pathway_info: dict {pathway_id: {'name': str, 'class': str}}
    """
    ko2path = {}
    pathway_info = {}
    for i, ko in enumerate(ko_list):
        try:
            if i > 0 and i % 10 == 0:
                time.sleep(0.2)
            entry = REST.kegg_get(ko).read()
            paths = []
            current_section = None
            for line in entry.split("\n"):
                if line.startswith("PATHWAY"):
                    parts = line.split()
                    if len(parts) > 1:
                        pid = parts[1]
                        paths.append(pid)
                        pname = " ".join(parts[2:]) if len(parts) > 2 else "Unknown"
                        if pid not in pathway_info:
                            pathway_info[pid] = {'name': pname, 'class': 'Unknown'}
                elif line.startswith("CLASS"):
                    current_section = "CLASS"
                elif current_section == "CLASS" and line.startswith(" "):
                    if paths:
                        pathway_info[paths[-1]]['class'] = line.strip()
            if paths:
                ko2path[ko] = paths
        except Exception as e:
            logging.warning(f"KO {ko} fetch error: {e}")
            continue
    return ko2path, pathway_info

def fisher_enrichment(ko2path, background_genes, pathway_info):
    study_counts = Counter([p for paths in ko2path.values() for p in paths])
    bg_ko2path, _ = get_ko_pathway_mapping(background_genes)
    bg_counts = Counter([p for paths in bg_ko2path.values() for p in paths])

    results = []
    total_study = len(ko2path)
    total_bg = len(background_genes)

    for pid, count in study_counts.items():
        bg_count = bg_counts.get(pid, 0)
        table = [[count, total_study - count],
                 [bg_count, total_bg - bg_count]]
        odds, pval = fisher_exact(table, alternative='greater')
        results.append({
            'Pathway': pid,
            'Pathway_name': pathway_info.get(pid, {}).get('name', 'Unknown'),
            'Pathway_class': pathway_info.get(pid, {}).get('class', 'Unknown'),
            'Study_count': count,
            'Background_count': bg_count,
            'Odds_ratio': odds,
            'P_value': pval
        })
    df = pd.DataFrame(results)
    if not df.empty:
        df['FDR'] = multipletests(df['P_value'], method='fdr_bh')[1]
        df.sort_values('FDR', inplace=True)
    return df

def plot_kegg_results(df, prefix, top_n=20):
    if df is None or df.empty:
        logging.warning(f"No KEGG results to plot for {prefix}")
        return
    top = df.head(top_n)

    # Barplot
    if 'FDR' in top.columns:
        plt.figure(figsize=(10,8))
        sns.barplot(x=-np.log10(top['FDR']), y=top['Pathway_name'], palette="viridis")
        plt.xlabel("-log10(FDR)")
        plt.ylabel("KEGG Pathway")
        plt.title(f"Top {top_n} KEGG Pathways ({prefix})")
        plt.tight_layout()
        plt.savefig(f"{prefix}_KEGG_barplot.png", dpi=300)
        plt.close()

    # Bubble plot
    if 'FDR' in top.columns and 'Odds_ratio' in top.columns:
        plt.figure(figsize=(10,8))
        sizes = top['Study_count'] * 20
        scatter = plt.scatter(top['Odds_ratio'], -np.log10(top['FDR']),
                              s=sizes, c=top['Study_count'], cmap='viridis', alpha=0.7)
        plt.colorbar(scatter, label='Gene count')
        plt.xlabel("Odds ratio")
        plt.ylabel("-log10(FDR)")
        plt.title(f"KEGG Bubble Plot ({prefix})")
        for i, row in top.iterrows():
            plt.text(row['Odds_ratio'], -np.log10(row['FDR']), row['Pathway'][:6])
        plt.tight_layout()
        plt.savefig(f"{prefix}_KEGG_bubbleplot.png", dpi=300)
        plt.close()

def kegg_enrichment(ko_list, background_genes, prefix="UP", top_n=20):
    """
    KEGG enrichment using Biopython REST + Fisher exact test.
    - ko_list: list of KOs to test
    - background_genes: list of KOs for background (required)
    - prefix: string prefix for output files
    - top_n: number of top pathways to plot
    Returns:
        df_summary: KEGG enrichment table (counts + FDR if background provided)
        df_extended: table mapping KOs to pathways
    """
    if not ko_list:
        logging.warning(f"No KOs provided for {prefix}")
        return None, None

    if not background_genes:
        raise ValueError("background_genes must be provided for KEGG enrichment.")

    # Map KOs to pathways
    ko2path, pathway_info = get_ko_pathway_mapping(ko_list)
    if not ko2path:
        logging.warning(f"No pathways found for {prefix} via KEGG REST")
        return None, None

    # Fisher enrichment
    df_summary = fisher_enrichment(ko2path, background_genes, pathway_info)

    # Save summary table
    df_summary.to_csv(f"{prefix}_KEGG_summary.tsv", sep="\t", index=False)
    plot_kegg_results(df_summary, prefix, top_n=top_n)

    # Extended KO -> pathway table
    extended_rows = []
    for ko, paths in ko2path.items():
        for p in paths:
            extended_rows.append({
                'KO': ko,
                'Pathway': p,
                'Pathway_name': pathway_info.get(p, {}).get('name','Unknown'),
                'Pathway_class': pathway_info.get(p, {}).get('class','Unknown')
            })
    df_extended = pd.DataFrame(extended_rows)
    df_extended.to_csv(f"{prefix}_KEGG_extended_table.tsv", sep="\t", index=False)

    return df_summary, df_extended

# ============================
# MAIN FUNCTION
# ============================
def main():
    args = parse_arguments()

    # Config logging
    configure_logging(args.logfile)

    # Validate inputs
    validate_inputs([args.ko_file, args.go_file, args.de_file, args.background, args.obo_file])

    # Load inputs
    ko_df, background, go_df, de_df = load_inputs(args)

    # Map gene names
    names_df = map_gene_names_to_new_df(ko_df, outdir=args.outdir, save_files=args.save_files)

    de_sets = split_up_down(
        de_df,
        fc_thresholds=args.fc_thresholds,
        alpha_values=args.alpha_values,
        fc_col='log2FoldChange',
        pval_col='padj'
    )

    for i, s in enumerate(de_sets):
        prefix = f"Set{i+1}_FC{str(s['fc']).replace('.', '_')}_alpha{str(s['alpha']).replace('.', '_')}"
        logging.info(f"Processing {prefix} ...")

        up_genes = s['up']
        down_genes = s['down']

        # Extract KOs and GOs
        files = extract_up_down_kos_gos(up_genes, down_genes, ko_df, go_df, outdir=args.outdir, save_files=args.save_files)

        # GO enrichment
        go_up, go_up_ext = go_enrichment(up_genes, go_df, de_df, background, args.obo_file, prefix=prefix + "_UP")
        go_down, go_down_ext = go_enrichment(down_genes, go_df, de_df, background, args.obo_file, prefix=prefix + "_DOWN")

        # KEGG enrichment
        kegg_up, kegg_up_ext = kegg_enrichment(list(files['up_KO']),
                                               list(files['up_KO']) + list(files['down_KO']),
                                               prefix=prefix + "_UP",
                                               top_n=args.top_n)
        kegg_down, kegg_down_ext = kegg_enrichment(list(files['down_KO']),
                                                   list(files['up_KO']) + list(files['down_KO']),
                                                   prefix=prefix + "_DOWN",
                                                   top_n=args.top_n)

    logging.info("Pipeline finished successfully.")

if __name__ == "__main__":
    main()
