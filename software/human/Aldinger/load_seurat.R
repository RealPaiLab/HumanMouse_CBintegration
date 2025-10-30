# ==============================================================================
# Loads the Aldinger RDS files
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/CBL_scRNAseq/software/")

library(dplyr)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
data_dir <- file.path("/isilon", root_dir, "data/human/Aldinger")
out_dir <- file.path("", root_dir, "results/human/Aldinger")
date_dir <- file.path(out_dir, "20220303")

# load the original Seurat object
orig_srat <- readRDS(file.path(data_dir, "seurat.rds"))

# load the Seurat object that Liam (Taylor lab) sent over
liam_srat <- readRDS(file.path(data_dir, "glutamatergic_dev_Liam.RDS"))

