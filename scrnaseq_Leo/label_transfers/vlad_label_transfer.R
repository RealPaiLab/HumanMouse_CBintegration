library(Seurat)
library(tidyverse)


dataset_srat <- readRDS(file.path("/data/llau/Vladoiu_2019/merged_seurat_RLonly.rds"))

# add cell names to mice cells
source("/u/llau/software/mb_scrnaseq/MB_scRNAseq/software/mouse/Vladoiu/add_annotations.R")
dataset_srat <- label_vladoiu_cells(dataset_srat)

# label transfer on unknown mice cells
dataset_srat <- SCTransform(dataset_srat, vars.to.regress = "CC.Difference")

# Subset with non-NA values 
dataset_srat.reference <- subset(dataset_srat, subset = mouse_cell_type != "NA - missing")

# Subset with NA values 
dataset_srat.query <- subset(dataset_srat, subset = mouse_cell_type == "NA - missing")

source("label_transfers/label_transfer.R")
label_list <- label_transfer(reference = dataset_srat.reference, query = dataset_srat.query, predict_column = "mouse_cell_type")

saveRDS(object = label_list[[1]], file = file.path("/data/llau/Vladoiu_2019/merged_seurat_RLonly_annotated.rds"))
write.csv(label_list[[2]], file = "/data/llau/Vladoiu_2019/merged_seurat_RLonly_annotated.csv", row.names = TRUE)