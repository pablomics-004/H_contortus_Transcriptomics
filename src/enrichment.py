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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
from math import isfinite
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
    de_df = pd.read_csv(args.de_file, sep="\t", header=0)
    de_df.columns = [c.strip() for c in de_df.columns]
    de_df = de_df.rename(columns={'Unnamed: 0': 'gene_id'})
    print(de_df.columns.tolist())
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
    
    # Si solo hay un threshold, usar alpha correspondiente
    if len(fc_thresholds_sorted) == 1:
        fc = fc_thresholds_sorted[0]
        alpha = alpha_values[0] if alpha_values else 0.05
        
        up = set(de_sub.loc[(de_sub[fc_col] >= fc) & (de_sub[pval_col] < alpha), 'gene_id'])
        down = set(de_sub.loc[(de_sub[fc_col] <= -fc) & (de_sub[pval_col] < alpha), 'gene_id'])
        
        logging.info(f"Single set: FC >= {fc}, alpha={alpha} -> Up: {len(up)}, Down: {len(down)}")
        results.append({'fc': fc, 'alpha': alpha, 'up': up, 'down': down})
    else:
        # Múltiples thresholds
        for i, (fc, alpha) in enumerate(zip(fc_thresholds_sorted, alpha_values)):
            if i == 0:
                # Primer set: genes entre 1er threshold y siguiente
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


def go_enrichment(genes_set: set, go_df: pd.DataFrame, de_df: pd.DataFrame, background: set, obo_file: str, prefix: str, outdir=None):
    """
    GO enrichment y tabla extendida. Devuelve (results_df, extended_df).
    - genes_set: set de gene_ids (strings)
    - go_df: DataFrame con columnas 'gene_id','GO'
    - de_df: DataFrame DE con columna 'gene_id' y métricas
    - background: iterable de gene_ids
    - obo_file: ruta al OBO
    - prefix: prefijo para archivos
    - outdir: directorio opcional donde guardar archivos
    """
    # preparar mapping gene -> GO (sin NaNs)
    gene2go = go_df.groupby('gene_id')['GO'].apply(lambda s: set([x for x in s if pd.notna(x)])).to_dict()

    # cargar ontología
    obodag = GODag(obo_file)

    # GOA expects population as list
    population = list(map(str, background))

    goea = GOEnrichmentStudy(population, gene2go, obodag,
                             propagate_counts=True, alpha=0.05, methods=['fdr_bh'])

    # run study (convert genes to str)
    study_genes = list(map(str, genes_set))
    results_all = goea.run_study(study_genes)

    rows = []
    extended = []

    for r in results_all:
        # some GOEnrichmentStudy results may not have p_fdr_bh attribute depending on version
        p_fdr = getattr(r, "p_fdr_bh", None)
        if p_fdr is None:
            # try alternative attribute names
            p_fdr = getattr(r, "p_fdr", None) or getattr(r, "pvalue_fdr_bh", None)

        if p_fdr is None:
            continue

        if p_fdr < 0.05:
            # obtener namespace seguro
            ns = getattr(r, "NS", None)
            if ns is None:
                ns = getattr(r, "namespace", None)
            # si viene como dict/OrderedDict, extraer campo razonable
            if isinstance(ns, (dict, list)):
                try:
                    if isinstance(ns, dict):
                        ns = ns.get("namespace") or ns.get("id") or ns.get("name")
                    else:
                        ns = str(ns)
                except Exception:
                    ns = str(ns)
            # mapear nombres largos a códigos BP/MF/CC
            ns_map_full = {
                'biological_process': 'BP',
                'molecular_function': 'MF',
                'cellular_component': 'CC',
                'Biological Process': 'BP',
                'Molecular Function': 'MF',
                'Cellular Component': 'CC'
            }
            ns_code = ns_map_full.get(str(ns), str(ns))

            # calcular ratios con seguridad
            ratio_in_study = np.nan
            ratio_in_pop = np.nan
            try:
                if getattr(r, 'ratio_in_study', None):
                    ratio_in_study = r.ratio_in_study[0] / r.ratio_in_study[1]
                if getattr(r, 'ratio_in_pop', None):
                    ratio_in_pop = r.ratio_in_pop[0] / r.ratio_in_pop[1]
            except Exception:
                pass

            rows.append({
                'GO': r.GO,
                'name': getattr(r, 'name', ''),
                'namespace': ns_code,
                'p_fdr': float(p_fdr),
                'ratio_in_study': ratio_in_study,
                'ratio_in_pop': ratio_in_pop,
                'n_genes': len(r.study_items) if getattr(r, 'study_items', None) else 0
            })

            # tabla extendida: listar genes del estudio dentro del término
            study_items = getattr(r, 'study_items', [])
            genes_in_term = [g for g in study_items if str(g) in set(map(str, genes_set))]
            for g in genes_in_term:
                # buscar fila en de_df si existe
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

    # guardar archivos
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

    # graficar si hay resultados
    if not results_df.empty:
        plot_go_enrichment(results_df, prefix, top_n=20, outdir=outdir)
    else:
        logging.warning(f"No significant GO terms found for {prefix}")

    return results_df, extended_df

def plot_go_enrichment(results: pd.DataFrame, prefix: str, top_n: int = 20, outdir=None):
    """
    Grafica principal: barh (-log10 FDR) + bubble plot (enrichment ratio)
    """
    if results is None or results.empty:
        print(f"No results for {prefix}")
        return

    # normalizar namespace a BP/MF/CC si es posible
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

    top = results.sort_values('p_fdr').head(top_n).copy()
    
    # DEBUG DETALLADO: Ver TODOS los datos
    print(f"\n=== GO RESULTS {prefix} ===")
    print(f"Total terms: {len(results)}, Top terms: {len(top)}")
    print(f"FDR range: {results['p_fdr'].min():.2e} to {results['p_fdr'].max():.2e}")
    print(f"Ratio_in_study range: {results['ratio_in_study'].min():.4f} to {results['ratio_in_study'].max():.4f}")
    print(f"Ratio_in_pop range: {results['ratio_in_pop'].min():.4f} to {results['ratio_in_pop'].max():.4f}")
    
    # Calcular enrichment ratio
    top['enrichment_ratio'] = top['ratio_in_study'] / top['ratio_in_pop']
    top['enrichment_ratio'] = top['enrichment_ratio'].replace([np.inf, -np.inf], np.nan)
    
    print(f"Enrichment ratio range: {top['enrichment_ratio'].min():.4f} to {top['enrichment_ratio'].max():.4f}")
    print(f"Gene counts range: {results['n_genes'].min()} to {results['n_genes'].max()}")
    
    # Mostrar los primeros 5 términos con sus ratios
    print("\nTop 5 terms with ratios:")
    for i, row in top.head().iterrows():
        print(f"  {row['name'][:50]}... | FDR: {row['p_fdr']:.2e} | Study: {row['ratio_in_study']:.4f} | Pop: {row['ratio_in_pop']:.4f} | Enrich: {row.get('enrichment_ratio', 'NaN'):.4f}")
    
    top['neg_log10_fdr'] = -np.log10(top['p_fdr'].replace(0, np.nextafter(0, 1)))

    # asignar colores por namespace
    ns_colors = {'BP': "#1d96ed", 'MF': "#f62020", 'CC': "#24e824"}
    top['color'] = top['namespace'].map(lambda x: ns_colors.get(x, '#7f7f7f'))

    # recortar nombres largos
    top['short_name'] = top['name'].apply(lambda s: (s[:80] + '...') if isinstance(s, str) and len(s) > 80 else s)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

    # --- GRÁFICA 1: Barras horizontales (-log10 FDR) ---
    if len(top) > 0:
        bars = ax1.barh(y=np.arange(len(top)), width=top['neg_log10_fdr'], color=top['color'])
        ax1.set_yticks(np.arange(len(top)))
        ax1.set_yticklabels(top['short_name'], fontsize=9)
        ax1.set_xlabel('-log10(FDR)')
        ax1.set_title(f"Top {len(top)} GO Terms — {prefix}", fontweight='bold')
        
        # AJUSTAR LÍMITES DEL EJE X para gráfica de barras
        max_fdr = top['neg_log10_fdr'].max()
        ax1.set_xlim(0, max_fdr * 1.15)
        
        print(f"Bar plot X limits: 0 to {max_fdr * 1.15:.2f}")
        
        # añadir etiquetas con n_genes y namespace
        for i, (_, row) in enumerate(top.iterrows()):
            txt = f"n={int(row.get('n_genes', 0))} | {row.get('namespace','')}"
            ax1.text(row['neg_log10_fdr'] + max_fdr * 0.01, i, txt, va='center', fontsize=8)
    else:
        ax1.text(0.5, 0.5, 'No significant terms', 
                transform=ax1.transAxes, ha='center', va='center', fontsize=12)
        ax1.set_title(f"GO Terms — {prefix}")

    # --- GRÁFICA 2: Bubble plot (Enrichment Ratio) ---
    # FILTRAR valores NaN para el bubble plot y RESETEAR ÍNDICES
    bubble_data = top.dropna(subset=['enrichment_ratio']).copy().reset_index(drop=True)
    
    if len(bubble_data) > 0:
        # CALCULAR TAMAÑOS desde bubble_data
        bubble_sizes = (bubble_data['n_genes'].fillna(1).astype(float)) * 25
        
        sc = ax2.scatter(bubble_data['enrichment_ratio'], np.arange(len(bubble_data)), 
                        s=bubble_sizes,
                        c=bubble_data['neg_log10_fdr'],
                        cmap='viridis', alpha=0.8, edgecolor='k')
        
        ax2.set_yticks(np.arange(len(bubble_data)))
        ax2.set_yticklabels(bubble_data['short_name'], fontsize=9)
        
        # AJUSTAR LÍMITES DEL EJE X para bubble plot
        min_ratio = bubble_data['enrichment_ratio'].min()
        max_ratio = bubble_data['enrichment_ratio'].max()
        
        # Asegurar márgenes adecuados
        x_margin = (max_ratio - min_ratio) * 0.1
        if x_margin == 0:  # Si todos los valores son iguales
            x_margin = max_ratio * 0.1 if max_ratio != 0 else 0.1
            
        ax2.set_xlim(min_ratio - x_margin, max_ratio + x_margin)
        
        # DEBUG: Ver los límites aplicados
        print(f"Bubble plot X limits: {min_ratio - x_margin:.3f} to {max_ratio + x_margin:.3f}")
        
        # Colorbar
        cbar = plt.colorbar(sc, ax=ax2)
        cbar.set_label('-log10(FDR)')
        
        # Línea de referencia
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

    # generar panel pequeño por namespace
    plot_go_by_namespace(results, prefix, top_n=10, outdir=outdir)

def plot_go_by_namespace(results: pd.DataFrame, prefix: str, top_n: int = 10, outdir=None):
    """
    Genera un PNG con 3 subplots (BP / MF / CC) mostrando -log10(FDR) para cada categoría.
    """
    if results is None or results.empty:
        return

    # normalizar namespace (misma lógica que antes)
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

    fig, axes = plt.subplots(1, 3, figsize=(22, 6), sharey=False)

    for ax, ns in zip(axes, namespaces):
        sub = results[results['namespace'] == ns].sort_values('p_fdr').head(top_n).copy()
        if sub is None or sub.empty:
            ax.set_title(f"{ns}: no terms")
            ax.axis('off')
            continue
        
        sub['neg_log10_fdr'] = -np.log10(sub['p_fdr'].replace(0, np.nextafter(0, 1)))
        sub['short_name'] = sub['name'].apply(lambda s: (s[:50] + '...') if isinstance(s, str) and len(s) > 50 else s)
        
        # CREAR GRÁFICA DE BARRAS
        bars = ax.barh(y=np.arange(len(sub)), width=sub['neg_log10_fdr'], color=colors.get(ns, '#7f7f7f'))
        
        # AJUSTAR LÍMITES DEL EJE X para cada subplot individualmente
        max_fdr_sub = sub['neg_log10_fdr'].max()
        ax.set_xlim(0, max_fdr_sub * 1.15)
        
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

def _chunked_iterable(iterable, size):
    it = list(iterable)
    for i in range(0, len(it), size):
        yield it[i:i+size]

def _parse_kegg_entry_text(entry_text):
    """
    Parse a single KEGG entry block text (from REST.kegg_get chunk) and return
    a tuple (entry_id, pathway_list, class_info_dict, name)
    This is a lightweight parser that reads PATHWAY and CLASS sections.
    """
    entry_id = None
    paths = []
    current_section = None
    last_path = None
    title = ""
    for line in entry_text.splitlines():
        if line.startswith("ENTRY"):
            parts = line.split()
            if len(parts) >= 2:
                entry_id = parts[1]
        elif line.startswith("NAME"):
            title = line.replace("NAME", "").strip()
        elif line.startswith("PATHWAY"):
            # PATHWAY   map00010  Glycolysis / Gluconeogenesis
            parts = line.split()
            if len(parts) >= 2:
                pid = parts[1]
                paths.append(pid)
                last_path = pid
        elif line.startswith("CLASS"):
            current_section = "CLASS"
        elif current_section == "CLASS" and line.startswith(" "):
            # class lines are indented; assign to last_path if exists
            if last_path:
                # return the class string (strip)
                # but we don't necessarily have a class mapping per pathway; handled outside
                pass
        else:
            # reset section if blank or different
            if line.strip() == "":
                current_section = None
    return entry_id, paths, title

def get_ko_pathway_mapping(ko_list, chunk_size=20, sleep=0.2, organism=None):
    """
    Map KO IDs to pathways using KEGG REST efficient multi-fetch.
    Params:
      - ko_list: iterable of KO identifiers (e.g. ['K00001','K00002'] or 'ko:K00001' variants)
      - chunk_size: number of KOs per REST.kegg_get call (tuneable)
      - sleep: pause between chunk requests
      - organism: optional KEGG organism code to *filter* pathways (e.g. 'cel')
    Returns:
      - ko2path: dict {ko_id_str: [pathway_id_str, ...]}
      - pathway_info: dict {pathway_id_str: {'name': str, 'class': str}}
    Notes:
      - This function normalizes KO ids (strip 'ko:' or 'K:' prefixes).
      - Filtering by organism is done by matching numeric pathway suffix (see code).
    """
    # normalize KO ids as strings
    kos = [str(k).strip() for k in ko_list if pd.notna(k)]
    kos = list(dict.fromkeys(kos))  # preserve order + unique

    # normalize forms like "ko:K00001" or "K00001" -> prefer "K00001"
    def _norm_ko(x):
        return re.sub(r'^(ko:|K:)?', '', x, flags=re.IGNORECASE).strip()

    kos_norm = list(map(_norm_ko, kos))

    ko2path = {}
    pathway_info = {}

    # Precompute valid pathway numeric suffixes for organism if requested
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

    cache = {}  # in-run cache: ko -> paths

    for chunk in _chunked_iterable(kos_norm, chunk_size):
        query = "+".join(chunk)
        try:
            text = REST.kegg_get(query).read()
        except Exception as e:
            logging.warning(f"KEGG REST multi-get failed for chunk starting with {chunk[0]}: {e}")
            # fallback to single-get per KO
            text = ""
            for k in chunk:
                try:
                    text += REST.kegg_get(k).read() + "\n\n"
                except Exception as e2:
                    logging.warning(f"KEGG REST single-get failed for {k}: {e2}")
                time.sleep(sleep)
        # split into entry blocks by 'ENTRY' lines: simpler to split on '\n\n'
        blocks = re.split(r'\n(?=ENTRY\s)', text)
        for block in blocks:
            if not block.strip():
                continue
            eid, paths, title = _parse_kegg_entry_text(block)
            if not eid:
                # try to infer eid from first word
                firstline = block.splitlines()[0]
                eid = firstline.split()[1] if len(firstline.split())>1 else None
            # choose a canonical ko key from the block or from chunk
            ko_key_candidates = []
            if eid:
                # normalize eid to Kxxxxx if it looks like K...
                ko_key_candidates.append(re.sub(r'^(ko:|K:)?', '', eid, flags=re.IGNORECASE))
            # also check lines for "ENTRY" content inside block
            # use first candidate as key
            if ko_key_candidates:
                ko_key = ko_key_candidates[0]
            else:
                # fallback: skip block
                continue

            # parse PATHWAY lines inside block more robustly
            block_paths = []
            for line in block.splitlines():
                if line.startswith("PATHWAY"):
                    parts = line.split()
                    if len(parts) >= 2:
                        pid = parts[1]
                        block_paths.append(pid)
                        if pid not in pathway_info:
                            pname = " ".join(parts[2:]) if len(parts) > 2 else "Unknown"
                            pathway_info[pid] = {'name': pname, 'class': 'Unknown'}
            # optional: parse CLASS lines (simple approach)
            # find "CLASS" then next indented line as class for last path if any
            cur_class = None
            for i, line in enumerate(block.splitlines()):
                if line.startswith("CLASS"):
                    cur_class = ""
                    # subsequent indented lines until blank or next section
                    j = i+1
                    while j < len(block.splitlines()) and block.splitlines()[j].startswith(" "):
                        cur_class += block.splitlines()[j].strip() + " "
                        j += 1
                    cur_class = cur_class.strip()
                    # assign to all block_paths as an approximation
                    for p in block_paths:
                        if p in pathway_info:
                            pathway_info[p]['class'] = cur_class

            # assign after optional organism filtering (we'll filter later)
            cache[ko_key] = block_paths

        time.sleep(sleep)

    # Build final ko2path using the cache, with optional organism filter
    for k, paths in cache.items():
        if not paths:
            continue
        kept = []
        if valid_path_digits is not None:
            for p in paths:
                m = re.search(r'(\d+)$', str(p))
                if m and m.group(1) in valid_path_digits:
                    kept.append(p)
        else:
            kept = paths[:]
        if kept:
            ko2path[k] = kept

    # pathway_info may contain paths not in ko2path; filter to keep only those present
    if pathway_info:
        allowed = set([p for paths in ko2path.values() for p in paths])
        pathway_info = {p: pathway_info[p] for p in pathway_info if p in allowed}

    return ko2path, pathway_info


def fisher_enrichment(ko2path_study, background_genes, pathway_info, organism=None, chunk_size=20, sleep=0.2):
    """
    Compute Fisher enrichment but fetch background mapping efficiently.
    - ko2path_study: dict of study KO -> [pathway,...]  (already computed)
    - background_genes: iterable of KOs (strings)
    - pathway_info: pathway_info dict (may be filtered by organism already)
    Returns df_summary (with FDR)
    """
    # If background mapping not already available, fetch for background but reuse existing pathway_info if possible
    bg_kos = [str(x).strip() for x in background_genes if pd.notna(x)]
    bg_kos_norm = list(dict.fromkeys([re.sub(r'^(ko:|K:)?', '', x, flags=re.IGNORECASE) for x in bg_kos]))

    # To avoid extra REST calls when we already have mapping for study KOs,
    # fetch mapping for background KOs that are not present in ko2path_study
    missing_bg = [k for k in bg_kos_norm if k not in ko2path_study]
    bg_ko2path = {}
    if missing_bg:
        bg_ko2path, _ = get_ko_pathway_mapping(missing_bg, chunk_size=chunk_size, sleep=sleep, organism=organism)
    # merge with existing (study) mappings if background includes study KOs
    # create combined bg mapping by using ko2path_study when available, else bg_ko2path
    combined_bg_map = {}
    for k in bg_kos_norm:
        if k in ko2path_study:
            combined_bg_map[k] = ko2path_study[k]
        elif k in bg_ko2path:
            combined_bg_map[k] = bg_ko2path[k]
        else:
            combined_bg_map[k] = []

    # counts
    study_counts = Counter([p for paths in ko2path_study.values() for p in paths])
    bg_counts = Counter([p for paths in combined_bg_map.values() for p in paths])

    results = []
    total_study = len(ko2path_study)
    total_bg = len([k for k in combined_bg_map.keys()])

    for pid, count in study_counts.items():
        bg_count = bg_counts.get(pid, 0)
        # contingency: study in pathway vs not, background in pathway vs not
        table = [[count, total_study - count],
                 [bg_count, total_bg - bg_count]]
        # ensure valid counts
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
    if not df.empty:
        df['FDR'] = multipletests(df['P_value'], method='fdr_bh')[1]
        df.sort_values('FDR', inplace=True)
    return df


def plot_kegg_results(df, prefix, top_n=20, outdir=None):
    """
    Improved KEGG plotting:
      - barplot: Study_count (x) vs Pathway_name (y); color = -log10(FDR)
      - bubble plot: x = Odds_ratio, y = Pathway_name, size = Study_count, color = -log10(FDR)
    """
    if df is None or df.empty:
        logging.warning(f"No KEGG results to plot for {prefix}")
        return

    df = df.copy()
    # ensure numeric
    df['Study_count'] = pd.to_numeric(df['Study_count'], errors='coerce').fillna(0).astype(int)
    df['FDR'] = pd.to_numeric(df.get('FDR', df.get('P_value')), errors='coerce').fillna(1.0)
    df['neg_log10_fdr'] = -np.log10(df['FDR'].replace(0, np.nextafter(0,1)))

    top = df.head(top_n).copy()
    top['short_name'] = top['Pathway_name'].apply(lambda s: (s[:80] + '...') if isinstance(s, str) and len(s) > 80 else s)
    top = top.sort_values('Study_count', ascending=True)  # so barh shows largest on top after invert

    # BAR PLOT: Study_count sized, color by significance
    fig, ax = plt.subplots(figsize=(10, max(6, 0.3 * len(top))))
    cmap = plt.get_cmap('viridis')
    colors = cmap((top['neg_log10_fdr'] - top['neg_log10_fdr'].min()) / (top['neg_log10_fdr'].max() - top['neg_log10_fdr'].min() + 1e-9))

    bars = ax.barh(top['short_name'], top['Study_count'], color=colors)
    ax.set_xlabel('Number of KOs (Study_count)')
    ax.set_title(f"Top {len(top)} KEGG Pathways ({prefix})")
    # add colorbar for -log10(FDR) by making a ScalarMappable
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=top['neg_log10_fdr'].min(), vmax=top['neg_log10_fdr'].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('-log10(FDR)')

    # annotate bars with FDR and Odds ratio
    for i, (_, row) in enumerate(top.iterrows()):
        ax.text(row['Study_count'] + max(1, int(0.01 * top['Study_count'].max())), i,
                f"FDR={row['FDR']:.2e}\nOR={row['Odds_ratio']:.2g}",
                va='center', fontsize=8)

    plt.tight_layout()
    outpath = f"{prefix}_KEGG_barplot.png"
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, outpath)
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()

    # BUBBLE PLOT: Odds_ratio vs -log10(FDR), size ~ Study_count
    bubble = top.copy().reset_index(drop=True)
    bubble_fig, bubble_ax = plt.subplots(figsize=(10, max(6, 0.3 * len(bubble))))
    sizes = (bubble['Study_count'].astype(float) + 1.0) * 30
    sc = bubble_ax.scatter(bubble['Odds_ratio'].replace([np.inf, -np.inf], np.nan),
                           bubble.index,
                           s=sizes,
                           c=bubble['neg_log10_fdr'],
                           cmap='viridis',
                           alpha=0.8,
                           edgecolor='k')
    bubble_ax.set_yticks(bubble.index)
    bubble_ax.set_yticklabels(bubble['short_name'], fontsize=9)
    bubble_ax.set_xlabel('Odds ratio')
    bubble_ax.set_title(f"KEGG Bubble Plot ({prefix})")
    bubble_ax.axvline(1.0, color='red', linestyle='--', alpha=0.6)
    cbar2 = plt.colorbar(sc, ax=bubble_ax)
    cbar2.set_label('-log10(FDR)')

    # annotate points with pathway id (short)
    for idx, row in bubble.iterrows():
        pid = str(row['Pathway'])[:8]
        bubble_ax.text(row['Odds_ratio'], idx, "  "+pid, va='center', fontsize=7)

    outbubble = f"{prefix}_KEGG_bubbleplot.png"
    if outdir:
        outbubble = os.path.join(outdir, outbubble)
    plt.tight_layout()
    plt.savefig(outbubble, dpi=300, bbox_inches='tight')
    plt.close()


def kegg_enrichment(ko_list, background_genes, prefix="UP", top_n=20, organism="cel", chunk_size=20, sleep=0.2, outdir=None):
    """
    Full KEGG enrichment pipeline optimized:
      - Fetch KO->pathway for union of study+background in chunks
      - Filter pathways by organism (if provided)
      - Compute Fisher enrichment (background-aware)
      - Save tables and produce plots (bar + bubble)
    Parameters:
      - ko_list: list of study KOs (strings)
      - background_genes: list of background KOs (strings)
      - organism: KEGG organism code for filtering (e.g. 'cel'); if None, no filtering
      - chunk_size: chunk size for multi-fetch
      - sleep: pause between chunk requests
      - outdir: optional output directory
    Returns:
      - df_summary, df_extended
    """
    # Normalize and unique
    study_kos = [re.sub(r'^(ko:|K:)?', '', str(x).strip(), flags=re.IGNORECASE) for x in ko_list if pd.notna(x)]
    bg_kos = [re.sub(r'^(ko:|K:)?', '', str(x).strip(), flags=re.IGNORECASE) for x in background_genes if pd.notna(x)]
    union_kos = list(dict.fromkeys(study_kos + bg_kos))

    if len(union_kos) == 0:
        logging.warning("No KOs provided for KEGG enrichment.")
        return None, None

    # Fetch mapping for union (one pass)
    ko2path_all, pathway_info_all = get_ko_pathway_mapping(union_kos, chunk_size=chunk_size, sleep=sleep, organism=organism)

    # split study mapping
    ko2path_study = {k: v for k, v in ko2path_all.items() if k in study_kos}
    # compute enrichment with efficient function
    df_summary = fisher_enrichment(ko2path_study, bg_kos, pathway_info_all, organism=organism, chunk_size=chunk_size, sleep=sleep)

    if df_summary is None or df_summary.empty:
        logging.warning(f"No KEGG enrichment results for {prefix}")
        return None, None

    # save summary and extended mapping
    fname_summary = f"{prefix}_KEGG_summary.tsv"
    fname_extended = f"{prefix}_KEGG_extended_table.tsv"
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        df_summary.to_csv(os.path.join(outdir, fname_summary), sep="\t", index=False)
    else:
        df_summary.to_csv(fname_summary, sep="\t", index=False)

    # extended table from ko2path_study
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

    # plotting
    plot_kegg_results(df_summary, prefix, top_n=top_n, outdir=outdir)

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

    # Map gene names (opcional, comentado por ahora)
    # names_df = map_gene_names_to_new_df(ko_df, outdir=args.outdir, save_files=args.save_files)

    de_sets = split_up_down(
        de_df,
        fc_thresholds=args.fc_thresholds,
        alpha_values=args.alpha_values,
        fc_col='log2FoldChange',
        pval_col='padj'
    )

    # Preparar background KOs una sola vez
    background_kos = (
        ko_df[ko_df["gene_id"].isin(background)]["KO"]
        .dropna()
        .unique()
        .tolist()
    )

    if len(background_kos) == 0:
        logging.error("No background genes have KO annotations — cannot run KEGG.")
        sys.exit(1)

    for i, s in enumerate(de_sets):
        up_genes = s['up']
        down_genes = s['down']
        
        # SKIP EMPTY SETS
        if len(up_genes) == 0 and len(down_genes) == 0:
            logging.info(f"Skipping Set{i+1} - no genes found")
            continue

        prefix = f"Set{i+1}_FC{str(s['fc']).replace('.', '_')}_alpha{str(s['alpha']).replace('.', '_')}"
        logging.info(f"Processing {prefix} ...")

        print("UP:", len(up_genes))
        print("DOWN:", len(down_genes))
        print("Background:", len(background))

        # Extract KOs and GOs
        files = extract_up_down_kos_gos(up_genes, down_genes, ko_df, go_df, outdir=args.outdir, save_files=args.save_files)

        # GO enrichment (solo si hay genes)
        if len(up_genes) > 0:
            go_up, go_up_ext = go_enrichment(up_genes, go_df, de_df, background, args.obo_file, prefix=prefix + "_UP", outdir=args.outdir)
        else:
            logging.info(f"No UP genes for GO enrichment in {prefix}")
            
        if len(down_genes) > 0:
            go_down, go_down_ext = go_enrichment(down_genes, go_df, de_df, background, args.obo_file, prefix=prefix + "_DOWN", outdir=args.outdir)
        else:
            logging.info(f"No DOWN genes for GO enrichment in {prefix}")

        # KEGG enrichment (solo si hay KOs)
        if len(files['up_KO']) > 0:
            kegg_up, kegg_up_ext = kegg_enrichment(
                ko_list=list(files['up_KO']),
                background_genes=background_kos,
                prefix=prefix + "_UP",
                top_n=args.top_n,
                outdir=args.outdir
            )
        else:
            logging.info(f"No UP KOs for KEGG enrichment in {prefix}")
            
        if len(files['down_KO']) > 0:
            kegg_down, kegg_down_ext = kegg_enrichment(
                ko_list=list(files['down_KO']),
                background_genes=background_kos,
                prefix=prefix + "_DOWN", 
                top_n=args.top_n,
                outdir=args.outdir
            )
        else:
            logging.info(f"No DOWN KOs for KEGG enrichment in {prefix}")

    logging.info("Pipeline finished successfully.")

if __name__ == "__main__":
    main()

