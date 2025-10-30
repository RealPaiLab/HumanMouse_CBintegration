library(biomaRt)
library(Seurat)

setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")
source("scrnaseq_Leo/utilities/get_orthologous_genes.R")

dataset_srat <- readRDS("/data/llau/Vladoiu_2019/merged_seurat_annotated.rds")
ortho_gene_directory <- "/data/llau/Vladoiu_2019"

get_orthologous_genes (
    dataset_srat,
    ortho_gene_directory
)