library(Seurat)
library(tidyverse)
library(qs)

setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")
source("scrnaseq_Leo/utilities/standardize_dataset.R")

# Sepp mouse full cerebellum
standardize_dataset(
  dataset_srat = readRDS("/.mounts/labs/pailab/private/llau/data/Sepp_2024/sepp2024_mouse_sct.rds"),
  dataset_author = "Sepp",
  dataset_name = "Sepp_full_cerebellum",
  species = "mouse",
  regress_var = "CC.Difference",
  cell_type_column = "cell_type",
  out_directory = "/.mounts/labs/pailab/private/llau/data/Sepp_2024"
)