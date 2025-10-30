# ==============================================================================
# Run pySCENIC on integrated human UBC subclusters.
# ==============================================================================

import os, subprocess, pickle, re, math
import operator as op
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from cytoolz import compose

from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.export import add_scenic_metadata
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
from IPython.display import HTML

# ------------------------------------------------------------------------------
# parse input arguments

parser = argparse.ArgumentParser()

parser.add_argument(
    # CSV of raw counts
    '--counts_csv',
    default=None,
    required=True
)
parser.add_argument(
    # CSV of metadata
    '--metadata_csv',
    default=None,
    required=True
)
parser.add_argument(
    # CSV of UMAP embeddings from Seurat object
    '--umap_csv',
    default=None,
    required=True
)
parser.add_argument(
    # column name with cluster information in the metadata table
    '--cluster_column',
    default=None,
    required=True,
)
parser.add_argument(
    '--num_threads',
    default=None,
    required=True,
    type=int
)
parser.add_argument(
    # prefix attached to all output files
    '--dataset_id',
    default=None,
    required=True
)
parser.add_argument(
    # output directory
    '--out_dir',
    default=None,
    required=True
)
parser.add_argument(
    # separate directory for figures; if `None`, defaults to output directory
    '--fig_dir',
    default=None,
    required=False
)

# args = parser.parse_args([
#     '--metadata_csv', '/.mounts/labs/pailab/private/projects/HumanMouseUBC/integrated_human_ubc/20250319/human_ubcs.metadata.csv',
#     '--cluster_column', 'SCT_snn_res.0.5',
#     '--num_threads', '4',
#     '--dataset_id', 'test',
#     '--out_dir', '/.mounts/labs/pailab/private/projects/HumanMouseUBC/integrated_human_ubc/20250319'
# ])
args = parser.parse_args()

# ------------------------------------------------------------------------------
# files, directories, and settings

# Set scanpy verbosity
sc.settings.verbosity = 3

# Set maximum number of jobs
sc.settings.njobs = 32

# counts, metadata, and UMAP embeddings from `seurat_to_scenic.R`
COUNTS_FNAME = args.counts_csv
METADATA_FNAME = args.metadata_csv
UMAP_EMBED_FNAME = args.umap_csv

# output files and results
RESULTS_DIR = args.out_dir
FIGURES_DIR = args.fig_dir if args.fig_dir is not None else args.out_dir
DATASET_ID = args.dataset_id
ANNDATA_FNAME = os.path.join(RESULTS_DIR, f'{DATASET_ID}.h5ad')
EXP_QC_FNAME = os.path.join(RESULTS_DIR, f'{DATASET_ID}.qc_counts.csv')
ADJACENCIES_FNAME = os.path.join(RESULTS_DIR, f'{DATASET_ID}.adjacencies.csv')
MOTIFS_ENR_FNAME = os.path.join(RESULTS_DIR, f'{DATASET_ID}.enr_motifs.csv')
REGULONS_DAT_FNAME = os.path.join(RESULTS_DIR, f'{DATASET_ID}.regulons.dat')
AUCELL_MTX_FNAME = os.path.join(RESULTS_DIR, f'{DATASET_ID}.auc_mtx.csv')
RSS_FNAME = os.path.join(RESULTS_DIR, f'{DATASET_ID}.rss.csv')

NUM_THREADS = args.num_threads

# make directories
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# set figure path
sc.settings.figdir = FIGURES_DIR

# pySCENIC files

# auxiliary cisTarget files
AUXILIARY_DIR = '/.mounts/labs/pailab/src/pySCENIC/cisTarget_databases'
RANKING_DBS_FNAMES = list(map(
    lambda fn: os.path.join(AUXILIARY_DIR, fn), 
    ['hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
     'hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather']
))
MOTIFS_ANNOT_FNAME = os.path.join(AUXILIARY_DIR, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')
TFS_FNAME = os.path.join(AUXILIARY_DIR, 'allTFs_hg38.txt')

# ------------------------------------------------------------------------------
# functions

# modified from the cancer data sets tutorial:
# http://htmlpreview.github.io/?https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/SCENIC%20Protocopyscenic.transform0study%20-%20Cancer%20data%20sets.html
def derive_regulons(motifs, db_fnames=RANKING_DBS_FNAMES):

    # if multi-level index, drop the top level
    if motifs.columns.nlevels == 2:
        motifs.columns = motifs.copy().columns.droplevel(0)

    db_names = tuple(map(lambda x: os.path.splitext(os.path.basename(x))[0], db_fnames))

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f

    # For the creation of regulons we only keep the 10-species databases and the
    # activating modules. We also remove the enriched motifs for the modules
    # that were created using the method 'weight>50.0%' (because these modules
    # are not part of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains(*db_names), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=bool)]

    # We build regulons only using enriched motifs with a NES of 3.0 or higher;
    # we take only directly annotated TFs or TF annotated for an orthologous
    # gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter(
        lambda r: len(r) >= 10,
        df2regulons(motifs[(motifs['NES'] >= 3.0) & (
            (motifs['Annotation'] == 'gene is directly annotated')
            | (motifs['Annotation'].str.startswith('gene is orthologous to')
               & motifs['Annotation'].str.endswith('which is directly annotated for motif'))
        )])
    ))

    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))

# function to calculate Z-score given the cell type annotations
def regulon_zscore(adata, cell_annot):

    # get the metadata of the cells
    df_metadata = adata.obs

    # get the column names of the regulon
    regulon_column_names = list(df_metadata.select_dtypes('number').columns)
    regulon_column_names = list(filter(lambda s: s.startswith('Regulon('), regulon_column_names))

    # subset cell type and regulon scores
    df_scores = df_metadata[regulon_column_names + [cell_annot]]

    # calculate Z-score [(sample - mean) / standard deviation
    sample_score = df_scores.groupby(by=cell_annot).mean()
    mean_score = df_metadata[regulon_column_names].mean()
    sd = df_metadata[regulon_column_names].std()
    df_zscore = ((sample_score - mean_score) / sd).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
    df_zscore['regulon'] = list(map(lambda s: s[8:-1], df_zscore.regulon))

    return df_zscore

# function to plot regulon activity
def plot_regulon_activity(adata, palette, color, color_map, **kwargs):
    """Plot expression of TFs and their regulons."""
    sc.set_figure_params(frameon=False, dpi=600, fontsize=8)
    
    # plot data
    fig = sc.pl.umap(adata, color=color, palette=palette, color_map=color_map, return_fig=True, **kwargs)
    
    return(fig)

# ------------------------------------------------------------------------------
# load data, convert to AnnData object

# load counts
df_counts = pd.read_csv(COUNTS_FNAME, index_col=0)

# load metadata
df_metadata = pd.read_csv(METADATA_FNAME, index_col=0)
df_metadata[args.cluster_column] = 'UBC_' + df_metadata[args.cluster_column].astype(str)
df_metadata.dataset_name = df_metadata.dataset_name.str.removesuffix('_full_cerebellum_human')
age_order = sorted(df_metadata.age.unique(), key = lambda x: (int(x.split()[0])))
df_metadata.age = df_metadata.age.astype('category').cat.reorder_categories(age_order)

# need to transpose counts so that rows are cells and columns are genes
adata = ad.AnnData(df_counts.T.sort_index(), dtype=np.float32)

# add metadata to AnnData object
adata.obs = df_metadata.sort_index()

sc.pp.filter_genes(adata, min_cells=5, inplace=True)

# save to output file
print(f'***Saving AnnData object to {ANNDATA_FNAME}***')
adata.write_h5ad(ANNDATA_FNAME)

# save gene expression matrix (only if running from scratch)
print(f'***Saving filtered gene expression matrix to {EXP_QC_FNAME}***')
adata.to_df().to_csv(EXP_QC_FNAME)

# ------------------------------------------------------------------------------
# run pySCENIC

# infer GRNs
grn_cmd = [
    'pyscenic', 'grn', EXP_QC_FNAME, TFS_FNAME,
    '-o', ADJACENCIES_FNAME,
    '--seed', '42',
    '--num_workers', f'{NUM_THREADS}'
]
subprocess.run(grn_cmd)

# run cisTarget to predict regulons
ctx_cmd = ['pyscenic', 'ctx', ADJACENCIES_FNAME] + RANKING_DBS_FNAMES \
    + ['--annotations_fname', MOTIFS_ANNOT_FNAME,
       '--expression_mtx_fname', EXP_QC_FNAME,
       '-o', MOTIFS_ENR_FNAME,
       '--num_workers', f'{NUM_THREADS}',
       '--mask_dropouts']
subprocess.run(ctx_cmd)
df_motifs = load_motifs(MOTIFS_ENR_FNAME)

# get regulons
regulons = derive_regulons(df_motifs)
# save (pickle) the regulons
print(f'***Saving regulons to {REGULONS_DAT_FNAME}***')
with open(REGULONS_DAT_FNAME, 'wb') as f:
    pickle.dump(regulons, f)

# ------------------------------------------------------------------------------
# run AUCell

# check number of genes to be used in AUCell
genes_detected_per_cell = np.sum(adata.X>0, axis=1)
print(f'***Total number of cells in the dataset: {len(genes_detected_per_cell)}***')
percentiles = pd.Series(genes_detected_per_cell).quantile(q=[.01, .05, .10, .50, 1])
print(f'***5% of cells express at least {percentiles[0.05]:.2f} genes***')

# plot distribution of the number of genes expressed in each cell
fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=600)
sns.histplot(genes_detected_per_cell, bins='fd')
for i,x in enumerate(percentiles):
    fig.gca().axvline(x=x, ymin=0,ymax=1, color='red')
    ax.text(x=x, y=ax.get_ylim()[1], s=f'{int(x)} ({percentiles.index.values[i]*100}%)', color='red', rotation=30, size='x-small', rotation_mode='anchor')
ax.set_xlabel('# of genes')
ax.set_ylabel('# of cells')
fig.tight_layout()
fig.savefig(os.path.join(FIGURES_DIR, 'histogram_num_genes_per_cell.png'))

# load expression matrix
exp_mtx = pd.read_csv(EXP_QC_FNAME, index_col=0)

# AUCell
auc_mtx = aucell(exp_mtx=exp_mtx, signatures=regulons, seed=42, num_workers=NUM_THREADS)

# save AUCell output
print(f'***Saving AUCell scores to {AUCELL_MTX_FNAME}***')
auc_mtx.to_csv(AUCELL_MTX_FNAME)

# ------------------------------------------------------------------------------
# UMAP projection and clustering

# add SCENIC metadata to `adata` object
print('***Adding SCENIC metadata to AnnData object***')
add_scenic_metadata(adata, auc_mtx, regulons)

# import and load UMAP embeddings from Seurat
print('***Loading UMAP embeddings and adding to AnnData object***')
umap_embeds = pd.read_csv(UMAP_EMBED_FNAME, index_col=0)
umap_embeds = umap_embeds.loc[adata.obs_names.to_list(), :]
adata.obsm['X_umap'] = umap_embeds.to_numpy(dtype='float32', copy=True)

print(f'***Saving AnnData object to {ANNDATA_FNAME}***')
adata.write_h5ad(ANNDATA_FNAME)

# plot UMAP
sc.set_figure_params(frameon=False, dpi=600, fontsize=8, figsize=(6, 6))
fig = sc.pl.umap(adata, color=[args.cluster_column, 'dataset_name', 'age'], ncols=2, return_fig=True)
fig.savefig(os.path.join(FIGURES_DIR, 'umap.png'), dpi=600)

# ------------------------------------------------------------------------------
# calculate z-score for each regulon in each cluster

print(f'***Calculating the z-score for each cluster in {args.cluster_column}***')
cluster_zscore = regulon_zscore(adata, args.cluster_column)
cluster_zscore.to_csv(os.path.join(RESULTS_DIR, 'regulon_zscore.csv'))

print('***Top 20 regulons:***')
print(cluster_zscore.sort_values('Z', ascending = False).head(20))

# format data for `sns.heatmap`
df_zscore_heatmap = pd.pivot_table(
    data=cluster_zscore,
    values='Z',
    index='regulon',
    columns=args.cluster_column
)

# make heatmap of z-scores
ht = df_zscore_heatmap.shape[0] * 0.2 + 2
wd = df_zscore_heatmap.shape[1] * 0.2 + 2
fig, ax = plt.subplots(figsize=(wd, ht))
ax = sns.heatmap(
    df_zscore_heatmap,
    annot=True,
    fmt='.1f',
    annot_kws={'size': 6},
    linewidths=0.7,
    linecolor='gray',
    cbar=True,
    cbar_kws={'shrink': 0.25},
    square=True,
    cmap='RdBu_r'
)
fig.tight_layout()

# save plot
fig_fname = os.path.join(FIGURES_DIR, 'regulon_zscore.png')
print(f'***Saving {fig_fname}***')
fig.savefig(fig_fname, dpi=600)

# ------------------------------------------------------------------------------
# get regulon specificity score for regulons in each cluster

rss = regulon_specificity_scores(auc_mtx, adata.obs[args.cluster_column])
print(f'***Saving regulon specificity scores to {RSS_FNAME}***')
rss.to_csv(RSS_FNAME)

sns.set_theme() # reset seaborn theme
sns.set_theme(style='whitegrid', font_scale=0.5)

# plot RSS for all clusters
ubc_clusters = adata.obs[args.cluster_column].unique().sort_values()
fig, axes = plt.subplots(nrows=math.ceil(len(ubc_clusters)/3), ncols=3, figsize=(6, 10))
for clust, ax in zip(ubc_clusters, axes.flat):
    plot_rss(rss, cell_type=clust, top_n=8, ax=ax)
    ax.set_xlabel('')
fig.tight_layout(h_pad=1, w_pad=1)

# save figure
fig.savefig(os.path.join(FIGURES_DIR, 'regulon_specificity_scores.png'), dpi=600)

# settings to plot top expressed TFs/regulons
COLORS = [color['color'] for color in mpl.rcParams['axes.prop_cycle']]
top_n = 8 # how many TFs to plot

# make and save plots of top expressed TFs/regulons
for clust in ubc_clusters:
    top_tfs = rss.loc[clust, :].nlargest(top_n).index.tolist()

    # expression of top TFs
    fig = plot_regulon_activity(
        adata,
        palette=COLORS,
        color=[args.cluster_column] + top_tfs,
        color_map='viridis',
        ncols=3,
        title=[clust] + top_tfs
    )
    fig_fname = os.path.join(FIGURES_DIR, re.sub('[/\s]', '_', clust) + '.TFs.png')
    print(f'***Saving {fig_fname}***')
    fig.savefig(fig_fname, dpi=600)

    # activity of top regulons
    top_regulons = list(map(lambda x: 'Regulon(' + x + ')', top_tfs))
    fig = plot_regulon_activity(
        adata,
        palette=COLORS,
        color=[args.cluster_column] + top_regulons,
        color_map='viridis',
        ncols=3,
        title=[clust] + top_regulons
    )

    # replace spaces and slashes for cell type filename
    fig_fname = os.path.join(FIGURES_DIR, re.sub('[/\s]', '_', clust) + '.regulons.png')
    print(f'***Saving {fig_fname}***')
    fig.savefig(fig_fname, dpi=600)