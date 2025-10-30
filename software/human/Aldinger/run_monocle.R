# ==============================================================================
# Use Monocle 3 to run trajectory inference on human RL
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(tidyverse)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # Hendrikse Seurat object
  "--liam_srat",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # integrated Seurat object
  "--integ_srat",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # number of cores to use
  "--num_cores",
  default = 32,
  required = FALSE
)

args <- parser$parse_args()

# for interactive testing only
# args <- parser$parse_args(c(
#   "--out_dir", "/CBL_scRNAseq/results/human/Aldinger/20221209/",
#   "--num_cores", "32"
# ))

if (is.null(args$out_dir)) {
  stop("Argument for `out_dir` is missing; please provide an output directory")
} else {
  out_dir <- args$out_dir
}
num_cores <- as.integer(args$num_cores)

# import functions
source("/CBL_scRNAseq/software/utilities/cell_labelling.R")

# ------------------------------------------------------------------------------
# import Seurat objects and copy cluster numbers from the integrated Seurat

liam_srat <- readRDS(args$liam_srat)
integ_srat <- readRDS(args$integ_srat)

# add cell type to integ_clusters
integ_srat <- label_integ_clusters(integ_srat, integ_method = "cca")

# copy the integrated cluster numbers to the human dataset, making sure that the
# cells have the correct cluster labelled
liam_srat$integ_clusters <- integ_srat@meta.data[row.names(liam_srat@meta.data), "integ_clusters"]

# ------------------------------------------------------------------------------
# trajectory inference analysis with Monocle 3

# convert to cell_data_set object (see
# https://github.com/satijalab/seurat-wrappers and
# https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html)
cds <- as.cell_data_set(liam_srat)

# normalize and calculate PCA (`preprocess_cds`)
message("Running PCA (`preprocess_cds`)")
cds <- preprocess_cds(
  cds,
  method = "PCA", # default is PCA
  num_dim = 100, # default is 50
  norm_method = "log" # default is "log"
)

# show proportion of variance explained
ggsave(
  "pc_variance_explained.png",
  plot = plot_pc_variance_explained(cds),
  path = out_dir,
  width = 8,
  height = 4,
  units = "in",
  dpi = 600
)

# dimensionality reduction to UMAP
message("Running UMAP (`reduce_dimension`)")
cds <- reduce_dimension(
  cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA"
)

# cluster cells with defaults
message("Grouping cells into clusters (`cluster_cells`)")
cds <- cluster_cells(cds, reduction_method = "UMAP")

# learn trajectory (with temporary workaround)
# names(cds@clusters[["UMAP"]]$partitions) <- names(cds@clusters[["UMAP"]]$clusters)
cds <- learn_graph(cds, use_partition = FALSE)

# plot cells by trajectory
plt <- plot_cells(
  cds,
  cell_size = 0.75,
  label_cell_groups = FALSE,
  label_principal_points = TRUE
)
ggsave(
  "trajectory_principal_points.png",
  plot = plt,
  path = out_dir,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

# order the cells by pseudotime, setting the root ("time 0") as the RL-VZ cells
cds <- order_cells(
  cds,
  root_cells = WhichCells(liam_srat, expression = new_cell_type == "RL-VZ")
)

# plot cells by pseudotime
plt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  cell_size = 0.75,
  label_cell_groups = FALSE,
)
ggsave(
  "trajectory_pseudotime.png",
  plot = plt,
  path = out_dir,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)

# plot cells by Liam's cluster labels
plt <- plot_cells(
  cds,
  color_cells_by = "new_cell_type",
  cell_size = 0.75,
  label_cell_groups = FALSE
)
ggsave(
  "trajectory_liam_cluster.png",
  plot = plt,
  path = out_dir,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)

# plot cells by integration cluster labels
plt <- plot_cells(
  cds,
  color_cells_by = "integ_clusters",
  cell_size = 0.75,
  label_cell_groups = FALSE
) + 
  scale_colour_manual(
    values = DiscretePalette(
      n = length(levels(colData(cds)$integ_clusters)),
      palette = "polychrome"
    )
  )
ggsave(
  "trajectory_integ_cluster.png",
  plot = plt,
  path = out_dir,
  width = 10,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# differential expression between the clusters

# run differential expression with `graph_test`
message("Identifying differentially expressed genes (`graph_test`)")
de_test_results <- graph_test(
  cds,
  neighbor_graph = "knn",
  cores = num_cores
)

# order by lowest q-value and largest Moran's I
de_test_results <- arrange(de_test_results, q_value, desc(abs(morans_I))) %>%
  # add row name as column
  rownames_to_column(var = "gene")

# write to file
message(sprintf(
  "Saving differential gene expression results to %s",
  file.path(out_dir, "de_genes.csv")
))
write_csv(de_test_results, file = file.path(out_dir, "de_genes.csv"))

# ------------------------------------------------------------------------------
# differential expression of homologous vs non-homologous UBCs

# # subset just our UBCs of interest
# # - cluster 5 is the conserved UBCs that clustered with mouse UBCs
# # - cluster 16 is the non-homologous UBCs that did not cluster with mouse UBCs
# ubc_cds <- cds[, colData(cds)$integ_clusters %in% c(5, 16)]
# 
# # run differential expression with `graph_test`
# message("Identifying differentially expressed genes (`graph_test`)")
# de_test_results <- graph_test(
#   ubc_cds,
#   neighbor_graph = "principal_graph",
#   alternative = "greater",
#   cores = num_cores
# )
# 
# # order by lowest q-value and largest Moran's I
# de_test_results <- arrange(de_test_results, q_value, desc(abs(morans_I))) %>% 
#   # add row name as column
#   rownames_to_column(var = "gene")
# 
# # write to file
# message(sprintf(
#   "Saving differential gene expression results to %s",
#   file.path(out_dir, "de_genes.csv")
# ))
# write_csv(de_test_results, file = "de_genes.csv", path = out_dir)

# ------------------------------------------------------------------------------
# find modules of co-regulated genes

sig_genes <- de_test_results[de_test_results$q_value < 0.05, "gene"] %>% 
  na.omit(.)

# cluster genes into modules
gene_module_df <- find_gene_modules(
  cds[sig_genes,],
  resolution = c(10^seq(-2, -6)),
  cores = num_cores,
  preprocess_method = "PCA"
)
write_csv(gene_module_df, file = file.path(out_dir, "gene_modules.csv"))

# get info about each cell group
cell_group_df <- tibble(
  cell = rownames(colData(cds)),
  integ_clusters = colData(cds)[, "integ_clusters"],
  cell_type = colData(cds)[, "new_cell_type"]
)
write_csv(cell_group_df, file = file.path(out_dir, "cell_groups.csv"))

# make heatmap of each cell group
walk(
  .x = seq(2, ncol(cell_group_df)),
  .f = function(x, cds, gene_group_df, cell_group_df, out_dir) {
    # subset cell grouping
    cell_group_df <- cell_group_df[, c(1, x)]
    
    # aggregate gene expression of modules for heatmap
    agg_gene_mat <- aggregate_gene_expression(
      cds,
      gene_group_df = gene_group_df,
      cell_group_df = cell_group_df
    )
    
    # make heatmap of gene clusters vs cell clusters
    hm <- pheatmap::pheatmap(
      mat = agg_gene_mat,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "column",
      clustering_method = "ward.D2",
      silent = TRUE
    )
    ggsave(
      filename = paste0("heatmap_", colnames(cell_group_df)[2], ".png"),
      plot = hm,
      path = out_dir,
      width = 0.25 * dim(agg_gene_mat)[2] + 1,
      height = 0.25 * dim(agg_gene_mat)[1] + 1,
      units = "in",
      dpi = 600
    )
  },
  cds = cds,
  gene_group_df = gene_module_df,
  cell_group_df = cell_group_df,
  out_dir = out_dir
)

# temporary workaround for plotting gene expression on UMAP (error likely
# results from converting a Seurat object to a CellDataSet object)
rowData(cds)$gene_short_name <- rownames(cds)

# plot expression of the modules
module_expr <- plot_cells(
  cds,
  reduction_method = "UMAP",
  genes = gene_module_df,
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE
)
ggsave(
  "module_expression.png",
  plot = module_expr,
  path = out_dir,
  width = ceiling(sqrt(length(levels(gene_module_df$module)))) * 4 + 1,
  height = ceiling(sqrt(length(levels(gene_module_df$module)))) * 3 + 1,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------

message(sprintf("Saving %s", file.path(out_dir, "monocle_cds.rds")))
saveRDS(cds, file = file.path(out_dir, "monocle_cds.rds"))

print(sessionInfo())
