library(Seurat)
library(tidyverse)
library(qs)

setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")
source("scrnaseq_Leo/utilities/standardize_dataset.R")

# Vladoiu full cerebellum
standardize_dataset(
  dataset_srat = readRDS("/.mounts/labs/pailab/private/llau/data/Vladoiu_2019/merged_seurat_annotated.rds"),
  dataset_author = "Vladoiu",
  dataset_name = "Vladoiu_full_cerebellum",
  species = "mouse",
  regress_var = "CC.Difference",
  cell_type_column = "mouse_cell_type",
  out_directory = "/.mounts/labs/pailab/private/llau/data/Vladoiu_2019"
)

# Vladoiu RL
standardize_dataset(
  dataset_srat = readRDS("/.mounts/labs/pailab/private/llau/data/Vladoiu_2019/merged_seurat_RLonly_annotated.rds"),
  dataset_author = "Vladoiu",
  dataset_name = "Vladoiu_RL",
  species = "mouse",
  regress_var = "CC.Difference",
  cell_type_column = "mouse_cell_type",
  out_directory = "/.mounts/labs/pailab/private/llau/data/Vladoiu_2019"
)