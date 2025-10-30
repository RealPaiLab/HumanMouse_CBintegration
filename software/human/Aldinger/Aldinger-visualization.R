# data from Aldinger et al. 2021 (PMID 34140698)
# https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html

# for use in the scRNAseq_docker container
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(ggplot2)
library(ggpubr)
library(Seurat)

# set paths for the expression matrix and metadata files
expr.path <- "/isilon/CBL_scRNAseq/data/human/Aldinger/exprMatrix.tsv.gz"
meta.path <- "/isilon/CBL_scRNAseq/data/human/Aldinger/meta.tsv"

# set path to save output figures
results.path <- "../results/20220125/"

# import the expression matrix and the meta data
exprMatrix <- read.csv(expr.path, sep = "\t", row.names = 1)
metadata <- read.csv(meta.path, sep = "\t", row.names = 1)

# create Seurat object
srat <- CreateSeuratObject(exprMatrix, project = "CBL-dev", meta.data = metadata)

# convert gene expression to Z-score (scale and centre the data)
srat <- ScaleData(srat, vars.to.regress = "CC.Difference", features = rownames(srat))

# run PCA for dimensionality reduction
all.genes <- row.names(srat)
srat <- RunPCA(srat, features = all.genes, npcs = 100)

# plot principal components
VizDimLoadings(srat, dims = 1:9) &
  theme(axis.text=element_text(size=6))
DimPlot(srat, reduction = "pca")
ElbowPlot(srat, ndims = 100)

# # function to cluster and make t-SNE plot
# make.tsne <- function(srat.obj, ndims, res) {
#   srat.obj <- FindNeighbors(srat.obj, dims = 1:ndims) %>% 
#     FindClusters(resolution = res)
#   
#   message("Generating t-SNE for ", ndims, " dimensions and ", res, " resolution.")
#   srat.obj <- RunTSNE(srat.obj, dims = 1:ndims)
#   tsne.plot <- DimPlot(srat.obj, reduction = "tsne", 
#                        repel = TRUE, label = TRUE)
#   
#   return(list(tsne.plot))
# }
# 
# # function to create multiple t-SNEs with different clustering parameters
# create.tsnes <- function(srat.obj, ndims, resolution) {
#   plots <- list()
#   
#   # loop through each of the parameters and append the plot to a list
#   for (dim in ndims) {
#     for (res in resolution) {
#       plot <- make.tsne(srat.obj, dim, res)
#       plots <- append(plots, plot)
#     }
#   }
#   
#   # return the list of plots
#   return(plots)
# }
# 
# # make the plots
# plots <- create.tsnes(srat, c(10, 25, 75), c(0.4, 0.8, 1.5))
# 
# # t-SNEs generated from automatic clustering and no annotation
# plot.noannot <- ggarrange(plotlist = plots, ncol = 3, nrow = 3, labels = "AUTO", legend = "none")
# ggsave(filename = "tsne-annot_none.pdf", plot = plot.noannot, width = 12, height = 12, units = "in")
# 
# # annotated t-SNE using the metadata cluster names
# plot.annot <- ggarrange(DimPlot(srat, reduction = "tsne", group.by = "Cluster", repel = TRUE, label = TRUE), 
#                         ncol = 1, nrow = 1, labels = "AUTO", legend = "right", common.legend = TRUE)
# ggsave(filename = "tsne-annot_paper.pdf", plot = plot.annot, width = 12, height = 6, units = "in")

# ------------------------------------------------------------------------------

# make t-SNE and UMAP plots using various clustering parameters then save as PDF

all.dims <- c(10, 25, 50, 75)
all.res <- c(0.1, 0.2, 0.4, 0.8, 1.2, 1.5)

for (dim in all.dims) {
  for (res in all.res) {
    srat.temp <- srat %>%
      FindNeighbors(dims = 1:dim) %>%
      FindClusters(resolution = res)
    
    message("Generating t-SNE for ", dim, " dimensions and ", res, " resolution.")
    srat.temp <- RunTSNE(srat.temp, dims = 1:dim, num_threads = 32)
    tsne.plot <- DimPlot(srat.temp, reduction = "tsne", repel = TRUE, label = TRUE) + 
      NoLegend()
    tsne.file <- paste0("tsne_dim", dim, "_res", res, ".pdf")
    ggsave(filename = tsne.file, plot = tsne.plot, path = results.path,
           height = 7, width = 7, units = "in")
    
    message("Generating UMAP for ", dim, " dimensions and ", res, " resolution.")
    srat.temp <- RunUMAP(srat.temp, dims = 1:dim)
    umap.plot <- DimPlot(srat.temp, reduction = "umap", repel = TRUE, label = TRUE) + 
      NoLegend()
    umap.file <- paste0("umap_dim", dim, "_res", res, ".pdf")
    ggsave(filename = umap.file, plot = umap.plot, path = results.path,
           height = 7, width = 7, units = "in")
  }
}

# ------------------------------------------------------------------------------

# annotate t-SNE and UMAP plots using the Aldinger annotations and clustering parameters

dim = 75

srat <- RunTSNE(srat, dims = 1:dim, num_threads = 32)
tsne.annot <- DimPlot(srat, reduction = "tsne", group.by = "figure_clusters", 
                      label = TRUE, repel = TRUE)
ggsave(filename = "tsne_annot.pdf", plot = tsne.annot, path = results.path, 
       height = 7, width = 10, units = "in")

srat <- RunUMAP(srat, dims = 1:dim)
umap.annot <- DimPlot(srat, reduction = "umap", group.by = "figure_clusters", 
                      label = TRUE, repel = TRUE)
ggsave(filename = "umap_annot.pdf", plot = umap.annot, path = results.path, 
       height = 7, width = 10, units = "in")

# ------------------------------------------------------------------------------

# plot gene expression on t-SNE and UMAP

genes <- c("LMX1A", 
           "EOMES", 
           "PAX6", 
           "WLS", 
           "MKI67", 
           "OTX2", 
           "RBFOX3", 
           "RELN")

for (gene in genes) {
  gene.plot <- FeaturePlot(srat, features = gene, reduction = "umap")
  ggsave(filename = paste0("umap_", gene, ".pdf"), plot = gene.plot, path = results.path, 
         height = 4, width = 5, units = "in")
}
