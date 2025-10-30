# ==============================================================================
# Prepare Luo 2022 dataset for integration. See QC in 2024-03-07 notes.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # subset cells for testing
  "--testing",
  default = NULL,
  type = "integer",
  required = FALSE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--out_dir", "/CBL_scRNAseq/results/human/Luo/20240315/",
    "--testing", 25
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# load Seurat object

con <- gzfile("/isilon/CBL_scRNAseq-archived/data/src/neurodev-genomics/scRNAseq/Luo_2022/Luo2022_humanFetalCerebellum_scRNAseq/GSM5952337_Fetal_cerebellum_final_cluster_1.rds.gz")
srat <- readRDS(gzcon(con))
close(con)

# subset cells for testing (so the script runs faster)
if (!is.null(arg_list$testing)) {
  message(sprintf("***Subsetting Seurat object to %s%%***", arg_list$testing))
  
  # get sample of cell IDs
  cell_sample <- sample(
    x = colnames(srat),
    size = ncol(srat)*(arg_list$testing*0.01),
    replace = FALSE
  )
  
  # subset Seurat object
  srat <- subset(srat, cells = cell_sample)
}

# save current cell types to column and change factor levels
srat$orig_cell_type <- Idents(srat)[rownames(srat[[]])]
new_levels <- c(
  "NSC",
  "TCP",
  "Dev.Purkinje",
  "Purkinje",
  "GABA interneuron",
  "GCP_P",
  "GCP_1",
  "GCP_2",
  "Granule_cell",
  "UBC_P",
  "Diff.UBC",
  "UBC",
  "Microglia",
  "Meninge",
  "Brainstem",
  "Monocyte",
  "Endothelial",
  "Glial progenitor",
  "Red_blood_cell",
  "Pericyte",
  "Vascular cell",
  "Ependymal cell",
  "Oligodendrocyte"
)
srat$orig_cell_type <- fct_relevel(srat$orig_cell_type, levels = new_levels)
levels(srat$orig_cell_type)

# ------------------------------------------------------------------------------
# calculate CC.Difference

srat <- CellCycleScoring(
  srat,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = FALSE
)

srat@meta.data <- mutate(
  srat@meta.data,
  CC.Difference = S.Score - G2M.Score
)

# ------------------------------------------------------------------------------
# run SCTransform

srat <- SCTransform(
  srat,
  variable.features.n = 5000,
  vars.to.regress = "CC.Difference"
  # return.only.var.genes = FALSE
)

message("\n***Free up RAM***\n")
gc()

# ------------------------------------------------------------------------------
# run PCA, UMAP, clustering

# save original PCA and UMAP reductions to new slot so they don't get overwritten
srat@reductions$pca_orig <- srat@reductions$pca
srat@reductions$umap_orig <- srat@reductions$umap

# re-run PCA/UMAP with SCT assay
srat <- RunPCA(srat, assay = "SCT", npcs = 100)
srat <- RunUMAP(srat, assay = "SCT", dims = 1:25)
srat <- FindNeighbors(srat, dims = 1:25)
srat <- FindClusters(srat, resolution = 0.8)

# ------------------------------------------------------------------------------
# save Seurat object

srat_out_path <- file.path(arg_list$out_dir, "luo2022_sct.rds")
message(sprintf("\n***Saving processed Seurat object to %s***\n", srat_out_path))
saveRDS(object = srat, file = srat_out_path)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
