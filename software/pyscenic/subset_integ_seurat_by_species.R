# ==============================================================================
# Subset human or mouse cells from a Seurat object after integration and perform
# dimensionality reduction.
# ==============================================================================

library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # path to Seurat object
  "--srat_qs",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # which species to subset (defaults to human)
  "--species",
  default = "human"
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  # >>> for testing >>>
  arg_list <- parser$parse_args(c(
    "--srat_qs", "/.mounts/labs/pailab/private/llau/results/integrated/20240524/25_pc_without_luo/25_pc_rl.qs",
    "--species", "mouse",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/pyscenic/20240704/"
  ))
  # <<< for testing <<<
} else {
  arg_list <- parser$parse_args()
}

if (!arg_list$species %in% c("human", "mouse")) {
  stop("`--species` must be either `human` or `mouse`")
}

message(sprintf("***Saving files to %s***", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# load Seurat object

srat <- qs::qread(arg_list$srat_qs)

# set default assay
DefaultAssay(srat) <- "integrated"

# set default Seurat clusters
srat$seurat_clusters <- srat$snn_res.0.4 %>%
  factor(x = ., levels = str_sort(unique(.), numeric = TRUE))

# save variable features
var_features <- VariableFeatures(srat, assay = "integrated")

# also need the UBC subclusters to copy them over
ubc_srat <- qs::qread("/.mounts/labs/pailab/private/llau/results/integrated/20240527/ubc_subset.qs")

# add UBC subclusters to the metadata
md <- merge(
  x = srat[[]],
  y = dplyr::select(ubc_srat[[]], snn_res.0.3) %>%
    dplyr::rename(ubc_subclusters = snn_res.0.3),
  by = "row.names",
  all.x = TRUE
) %>%
  mutate(
    # combine UBC subclusters with other clusters
    final_clusters = case_when(
      !is.na(ubc_subclusters) ~ paste0("UBC_", ubc_subclusters),
      .default = seurat_clusters
    ),
    # set factor levels for the new column
    final_clusters = fct_relevel(
      final_clusters,
      str_sort(names(table(final_clusters)), numeric = TRUE)
    )
  ) %>%
  column_to_rownames(var = "Row.names")

# add cell type annotations to the metadata
cluster_annot <- read_csv("/.mounts/labs/pailab/private/llau/results/integrated/20240618/clusters_cell_types.csv") %>%
  mutate(
    final_clusters = str_remove(Cluster, "^Cluster_"),
    .before = 1
  ) %>%
  dplyr::rename(cell_type_annot = broad_w_ubc_subtypes) %>%
  dplyr::select(final_clusters, cell_type_annot)
md <- merge(
  x = rownames_to_column(md, var = "cell_id"),
  y = cluster_annot,
  by = "final_clusters",
  all.x = TRUE
) %>%
  column_to_rownames(var = "cell_id")

# keep same order as before and replace the metadata
md <- md[rownames(srat[[]]), ]
srat@meta.data <- md

# ------------------------------------------------------------------------------
# subset human cells

message(sprintf(
  "***Subsetting %s cells from the integrated Seurat object***",
  arg_list$species
))
srat <- subset(
  srat,
  subset = species == arg_list$species
)

# ------------------------------------------------------------------------------
# re-run dimensionality reduction

# PCA
message(sprintf(
  "***Running PCA on the subsetted %s cells***",
  arg_list$species
))
srat <- RunPCA(srat, npcs = 100)

# UMAP
message(sprintf(
  "***Running UMAP on the subsetted %s cells using the top 25 PCs***",
  arg_list$species
))
srat <- RunUMAP(srat, dims = 1:25)

# ------------------------------------------------------------------------------
# save object

# use RDS file cuz downstream script doesn't take qs::qsave files
fname <- paste0("integ_", arg_list$species, "_srat.rds")
message(sprintf(
  "***Saving the subsetted %s cells as %s***",
  arg_list$species,
  fname
))
saveRDS(object = srat, file = file.path(arg_list$out_dir, fname))

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
