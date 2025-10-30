library(Seurat)
library(tidyverse)


dataset_srat <- readRDS(file.path("/data/llau/Cao_2020/Cao2020_10k_subset_cerebellum_seurat.rds"))

dataset_srat <- SCTransform(dataset_srat)

dataset_srat@meta.data$Main_cluster_name[is.na(dataset_srat@meta.data$Main_cluster_name)] <- "Unknown"

# Subset with non-NA values 
dataset_srat.reference <- subset(dataset_srat, subset = Main_cluster_name != "Unknown")

# Subset with NA values 
dataset_srat.query <- subset(dataset_srat, subset = Main_cluster_name == "Unknown")

source("label_transfers/label_transfer.R")
label_list <- label_transfer(reference = dataset_srat.reference, query = dataset_srat.query, predict_column = "Main_cluster_name")

saveRDS(object = label_list[[1]], file = file.path("/data/llau/Cao_2020/Cao2020_10k_subset_cerebellum_seurat_annotated.rds"))
write.csv(label_list[[2]], file = "/data/llau/Cao_2020/Cao2020_10k_subset_cerebellum_seurat_annotated.csv", row.names = TRUE)