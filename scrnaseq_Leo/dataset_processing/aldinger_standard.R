library(Seurat)
library(tidyverse)
library(qs)

setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")
source("scrnaseq_Leo/utilities/standardize_dataset.R")

# Aldinger full cerebellum
standardize_dataset(
  dataset_srat = readRDS("/.mounts/labs/pailab/private/llau/data/Aldinger_2021/seurat.rds"),
  dataset_author = "Aldinger",
  dataset_name = "Aldinger_full_cerebellum",
  species = "human",
  regress_var = "CC.Difference",
  cell_type_column = "fig_cell_type",
  out_directory = "/.mounts/labs/pailab/private/llau/data/Aldinger_2021"
)

# Aldinger RL
standardize_dataset(
  dataset_srat = readRDS("/.mounts/labs/pailab/private/llau/data/Aldinger_2021/glutamatergic_dev_Liam.RDS"),
  dataset_author = "Aldinger",
  dataset_name = "Aldinger_RL",
  species = "human",
  regress_var = "CC.Difference",
  cell_type_column = "new_cell_type",
  out_directory = "/.mounts/labs/pailab/private/llau/data/Aldinger_2021"
)