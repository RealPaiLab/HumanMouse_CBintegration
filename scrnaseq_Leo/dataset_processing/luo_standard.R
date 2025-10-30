library(Seurat)
library(tidyverse)
library(qs)

setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")
source("scrnaseq_Leo/utilities/standardize_dataset.R")

# Luo full cerebellum
standardize_dataset(
  dataset_srat = readRDS("/.mounts/labs/pailab/private/llau/data/Luo_2022/luo2022_sct.rds"),
  dataset_author = "Luo",
  dataset_name = "Luo_full_cerebellum",
  species = "human",
  regress_var = "CC.Difference",
  cell_type_column = "orig_cell_type",
  out_directory = "/.mounts/labs/pailab/private/llau/data/Luo_2022"
)