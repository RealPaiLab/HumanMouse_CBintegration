# see Aldinger et al (PMID 34140698)

setwd("/data/scRNAseq/CBL-dev/")

library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratWrappers)
library(monocle3)

# import the expression matrix and the meta data
exprMatrix <- read.csv("exprMatrix.tsv.gz", sep = "\t", row.names = 1)
metadata <- read.csv("meta.tsv", sep = "\t", row.names = 1) %>% 
  mutate(Cluster = as.factor(Cluster), 
         figure_clusters = as.factor(figure_clusters), 
         age = as.factor(age))

# create Seurat object
srat <- CreateSeuratObject(exprMatrix, project = "CBL-dev", meta.data = metadata)

# convert gene expression to Z-score (scale and centre the data)
srat <- ScaleData(srat, vars.to.regress = "CC.Difference", features = rownames(srat))


# ------------------------------------------------------------------------------
# make trajectory of rhombic lip cells

# subset the cells that originate from the RL
rl_cells = c("02-RL", "03-GCP", "04-GN", "05-eCN/UBC")
rl.srat <- subset(srat, subset = Cluster %in% rl_cells)

# run principal components analysis and show elbow plot
all.genes <- row.names(rl.srat)
rl.srat <- RunPCA(rl.srat, features = all.genes, npcs = 75)
ElbowPlot(rl.srat, ndims = 75)

# run t-SNE and UMAP
ndims <- 20 # determined from ElbowPlot
# rl.srat <- RunTSNE(rl.srat, dims = 1:ndims)
rl.srat <- RunUMAP(rl.srat, dims = 1:ndims)

# clustering with Seurat
rl.srat <- FindNeighbors(rl.srat, dims = 1:ndims)
rl.srat <- FindClusters(rl.srat)
DimPlot(rl.srat)
DimPlot(rl.srat, group.by = "figure_clusters")

# create CellDataSet object from Seurat object
rl.cds <- as.cell_data_set(rl.srat)

# cluster the CellDataSet (Monocle3)
rl.cds <- preprocess_cds(rl.cds, method = "PCA", num_dim = ndims)
rl.cds <- reduce_dimension(rl.cds, reduction_method = "UMAP", preprocess_method = "PCA")
rl.cds <- cluster_cells(rl.cds, reduction_method = "UMAP")
plot_cells(rl.cds, reduction_method = "UMAP", show_trajectory_graph = FALSE, 
           color_cells_by = "figure_clusters", label_groups_by_cluster = FALSE, 
           cell_size = 0.5, group_label_size = 4)

# learn the trajectory graph
rl.cds <- learn_graph(rl.cds)
plot_cells(rl.cds, reduction_method = "UMAP", color_cells_by = "figure_clusters", 
           label_groups_by_cluster = TRUE, cell_size = 0.5, group_label_size = 4)

rl.cds <- order_cells(rl.cds)
plot_cells(rl.cds, reduction_method = "UMAP", color_cells_by = "pseudotime", 
           cell_size = 0.5, group_label_size = 4)
