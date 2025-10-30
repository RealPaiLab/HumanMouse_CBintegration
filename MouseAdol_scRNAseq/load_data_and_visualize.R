# Data can be found in "/.mounts/labs/pailab/private/icheong/MouseAdol_scRNAseq"
# PMID: 31519873 (Bhattacherjee et al. 2019 Nat Commun)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(ggplot2)
library(Seurat)

# set input and output directories
data_dir <- "/isilon/MouseAdol_scRNAseq"
out_dir <- "/CBL_scRNAseq/MouseAdol_scRNAseq/results"
date_dir <- file.path(out_dir, format(Sys.Date(), "%Y%m%d"))

# ------------------------------------------------------------------------------
# load matrix/metadata and make Seurat object

counts <- read.csv(file.path(data_dir, "GSE124952_expression_matrix.csv"), 
                   row.names = 1)
metadata <- read.csv(file.path(data_dir, "GSE124952_meta_data.csv"),
                     row.names = 1)

srat <- CreateSeuratObject(
  counts = counts,
  project = "Mouse_Adol",
  meta.data = metadata
)

# ------------------------------------------------------------------------------
# integrate the adult (saline controls) and P21 data with CCA

# make a list of two Seurat objects: Adult (saline controls) and P21
srat_age <- SplitObject(srat, split.by = "DevStage")
srat_treat <- SplitObject(srat_age$Adult, split.by = "treatment")
srat_list <- list(Adult = srat_treat$Saline, P21 = srat_age$P21)

# normalize data and find variable features (default parameters are explicitly listed)
srat_list <- lapply(
  X = srat_list,
  FUN = function(X) {
    X <- NormalizeData(X, normalization.method = "LogNormalize", scale.factor = 10000)
    X <- FindVariableFeatures(X, selection.method = "vst", nfeatures = 10000)
  }
)

# select common features for integration
features <- SelectIntegrationFeatures(object.list = srat_list, nfeatures = 10000)

# identify integration anchors
anchors <- FindIntegrationAnchors(srat_list, anchor.features = features)

# integration
integ <- IntegrateData(anchors, normalization.method = "LogNormalize")

# ------------------------------------------------------------------------------
# visualize integrated data

# scale all genes (paper does not mention any variables to regress out)
integ <- ScaleData(integ, features = rownames(integ[["RNA"]]), vars.to.regress = NULL)

# run PCA
integ <- RunPCA(integ, npcs = 50)

# use JackStraw to identify number of components for clustering
# paper used p < 1E-3 = 0.001
# integ <- JackStraw(integ, dims = 50, num.replicate = 100) %>% 
#   ScoreJackStraw(dims = 1:50)
# JackStrawPlot(integ, dims = 1:50)

ElbowPlot(integ, ndims = 50)
ggsave("elbow.pdf", path = date_dir, height = 6, width = 8, units = "in")

# use top 20 PCs for downstream analysis
ndim <- 20

# run t-SNE
integ <- RunTSNE(integ, dims = 1:ndim, num_threads = 32)

# run UMAP
integ <- RunUMAP(integ, dims = 1:ndim)

# save plots
for (reduction in c("tsne", "umap")) {
  for (group in c("CellType", "DevStage")) {
    ggsave(
      paste0(reduction, "_", group, ".pdf"),
      plot = DimPlot(integ, reduction = reduction, group.by = group),
      path = date_dir,
      width = 6,
      height = 5,
      units = "in"
    )
  }
}

# save integrated Seurat object
saveRDS(integ, file = file.path(dirname(out_dir), "untreated_srat.rds"))
