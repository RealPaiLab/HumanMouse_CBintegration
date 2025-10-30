# ==============================================================================
# Save Seurat raw counts and metadata as CSV for use in SCENIC.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # path to Seurat object
  "--srat_rds",
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
  # counts file
  "--out_counts",
  default = NULL,
  required = TRUE
)

parser$add_argument(
  # metadata file
  "--out_metadata",
  default = NULL,
  required = TRUE
)

parser$add_argument(
  # metadata file
  "--out_umap",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  args <- parser$parse_args(c(
    "--srat_rds", "/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS",
    "--out_dir", getwd(),
    "--out_counts", "glutamatergic_dev_Liam.rna_counts.csv",
    "--out_metadata", "glutamatergic_dev_Liam.metadata.csv",
    "--out_umap", "glutamatergic_dev_Liam.umap_embeddings.csv"
  ))
} else {
  args <- parser$parse_args()
}

# ------------------------------------------------------------------------------
# import data

message(sprintf("Reading in from: %s", args$srat_rds))
srat <- readRDS(args$srat_rds)

# message("Reading in Seurat object from: /CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds")
# integ_srat <- readRDS("/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds")

# ------------------------------------------------------------------------------
# import functions

# source("/CBL_scRNAseq/software/utilities/cell_labelling.R")

# ------------------------------------------------------------------------------
# save counts

counts_path <- file.path(args$out_dir, args$out_counts)
slot <- "counts"
counts <- GetAssayData(srat[["RNA"]], slot = slot)

message(sprintf("Saving counts to %s", counts_path))
write.csv(as.data.frame(counts), file = counts_path)

# ------------------------------------------------------------------------------
# save metadata

metadata_path <- file.path(args$out_dir, args$out_metadata)

# copy homologous/non-homologous UBC labels from the integrated dataset to the human dataset
# integ_srat <- label_integ_clusters(
#   integ_srat,
#   col_name = "integ_clusters",
#   integ_method = "cca",
#   datasets = "rl"
# )
# srat$integ_clusters <- integ_srat$integ_clusters[rownames(srat@meta.data)]

# srat[[]] retrieves metadata, same as srat@meta.data
message(sprintf("Saving metadata to %s", metadata_path))
write.csv(srat[[]], file = metadata_path)

# ------------------------------------------------------------------------------
# save UMAP embeddings

umap_path <- file.path(args$out_dir, args$out_umap)
umap_embeds <- srat[["umap"]]@cell.embeddings %>% as.data.frame()

message(sprintf("Saving UMAP embeddings to %s", umap_path))
write.csv(umap_embeds, file = umap_path)

# ------------------------------------------------------------------------------

message("\nSESSION INFO\n")
print(sessionInfo())
