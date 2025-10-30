# ==============================================================================
# Isolate the human UBCs then reprocess and recluster for downstream
# differential gene expression analysis. The UBCs are identified based on the
# integrated clustering of the human RL with the mouse RL.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(Seurat)

# set output directory
date_dir <- "/CBL_scRNAseq/results/human/Aldinger/20220715"

# ------------------------------------------------------------------------------
# get the UBC cell names from the integrated Seurat object

integ_srat <- readRDS("/CBL_scRNAseq/results/integrated/vladoiu_liam_RL.rds")

ubc_barcodes <-
  integ_srat@meta.data[integ_srat$species == "human" &
                         integ_srat$seurat_clusters %in% c(8, 17, 18), ] %>% 
  rownames(.)

# ------------------------------------------------------------------------------
# subset UBC cells from Aldinger data

aldin_srat <- readRDS("/isilon/CBL_scRNAseq/data/human/Aldinger/glutamatergic_dev_Liam.RDS")
aldin_srat$is_ubc <- rownames(aldin_srat@meta.data) %in% ubc_barcodes

# change active assay to RNA so that SCT can be deleted
aldin_srat@active.assay <- "RNA"

ubc_srat <- subset(
  x = aldin_srat,
  subset = is_ubc == TRUE
) %>% 
  DietSeurat(
    .,
    assays = "RNA"
  )

# ------------------------------------------------------------------------------
# copy cluster idents over from the integrated metadata to the Aldinger metadata

ubc_srat$integ_clusters <-
  integ_srat$seurat_clusters[rownames(integ_srat@meta.data) %in% ubc_barcodes] %>% 
  factor(.)

# ------------------------------------------------------------------------------
# reprocess UBCs by themselves

source("/CBL_scRNAseq/software/utilities/process_raw_seurat.R")

ubc_srat <- process_seurat(
  ubc_srat,
  out_path = date_dir,
  ndims = 30, # Quang said they use 30 as default
  RunTSNE_args = list(num_threads = 32)
)

saveRDS(object = ubc_srat, file = file.path(dirname(date_dir), "UBC_seurat.rds"))
