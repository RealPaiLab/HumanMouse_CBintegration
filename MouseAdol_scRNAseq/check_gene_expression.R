# Data can be found in "/.mounts/labs/pailab/private/icheong/MouseAdol_scRNAseq"
# PMID: 31519873 (Bhattacherjee et al. 2019 Nat Commun)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(Seurat)

# set input and output directories
data_dir <- "/isilon/MouseAdol_scRNAseq"
out_dir <- "/CBL_scRNAseq/MouseAdol_scRNAseq/results"
date_dir <- file.path(out_dir, format(Sys.Date(), "%Y%m%d"))

# get function to integrate Seurat object
source("../software/utilities/process_raw_seurat.R")

# ------------------------------------------------------------------------------
# load Seurat objects

srat <- readRDS(file.path(dirname(out_dir), "untreated_srat.rds"))
srat <- SetIdent(srat, value = "CellType")
srat$DevStage <- factor(srat$DevStage, levels = c("P21", "Adult"))

aldinger_srat <- readRDS("/isilon/CBL_scRNAseq/data/human/Aldinger/seurat.rds")

# ------------------------------------------------------------------------------
# get differentially expressed genes between timepoints within each cell type
# e.g. astrocytes P21 vs. adult, oligocytes P21 vs. adult, etc.

cell_types <- c("Astro", "Oligo")
DefaultAssay(srat) <- "RNA"

for (type in cell_types) {
  # make directory for cell type
  dir.create(path = file.path(date_dir, type), showWarnings = FALSE)
  
  # isolate cell type and split by developmental stage for re-integration
  srat_list <- subset(srat, subset = CellType == type) %>%
    DietSeurat(., assays = "RNA") %>% 
    SplitObject(., split.by = "DevStage")
  
  # integrate the two timepoints (custom function in `utilities` directory)
  cell_type_srat <- integrate_seurat(srat_list, nfeatures = 2000)
  
  # run PCA, UMAP
  ndims = 30
  cell_type_srat <- RunPCA(cell_type_srat, npcs = ndims)
  cell_type_srat <- RunUMAP(cell_type_srat, dims = 1:ndims)
  cell_type_srat <- FindNeighbors(cell_type_srat, reduction = "pca", dims = 1:ndims)
  cell_type_srat <- FindClusters(cell_type_srat, resolution = 0.5)
  
  # save UMAP
  umap_clusters <- DimPlot(cell_type_srat, reduction = "umap", label = TRUE, repel = TRUE) + 
    NoLegend()
  umap_age <- DimPlot(cell_type_srat, reduction = "umap", group.by = "DevStage", label = TRUE, repel = TRUE)
  ggsave(
    filename = paste0("umap_", type, ".png"),
    plot = umap_age + umap_clusters,
    path = file.path(date_dir, type),
    width = 8,
    height = 4,
    units = "in",
    dpi = 300
  )
  
  DefaultAssay(cell_type_srat) <- "RNA"
  
  # differential gene expression
  de_genes <- FindMarkers(
    cell_type_srat,
    # assay = "RNA",
    ident.1 = "Adult",
    ident.2 = "P21",
    group.by = "DevStage",
    test.use = "wilcox"
  ) %>% 
    mutate(gene = rownames(.), .before = 1)
  
  # save differential gene expression results
  write_csv(
    de_genes,
    file = file.path(date_dir, type, paste0("de_genes_", type, ".csv")),
    col_names = TRUE
  )
  
  # get differential expression of Fbln1 and Atxn10 only
  genes_of_interest <- FindMarkers(
    cell_type_srat,
    # assay = "RNA",
    ident.1 = "Adult",
    ident.2 = "P21",
    group.by = "DevStage",
    features = c("Atxn10", "Fbln1"),
    logfc.threshold = 0,
    min.pct = 0,
    test.use = "wilcox"
  ) %>% 
    mutate(gene = rownames(.), .before = 1)
  
  # save results
  write_csv(
    genes_of_interest,
    file = file.path(date_dir, type, "de_genes_of_interest.csv"),
    col_names = TRUE
  )
}

# ------------------------------------------------------------------------------
# average expression of Fbln1 and Atxn10 in P21 vs. Adult with Wilcoxon test

de_genes_of_interest <- FindMarkers(
  srat,
  assay = "RNA",
  ident.1 = "Adult",
  ident.2 = "P21",
  group.by = "DevStage",
  features = c("Atxn10", "Fbln1"),
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

# save these genes
write.csv(de_genes_of_interest, file = file.path(date_dir, "de_genes_of_interest.csv"), quote = FALSE)

# get average expression of these genes
avg_exp <- AverageExpression(srat, assays = "RNA", features = c("Atxn10", "Fbln1"), group.by = "DevStage")[[1]]
write.csv(avg_exp, file = file.path(date_dir, "avg_exp_genes_of_interest.csv"), quote = FALSE)

# ------------------------------------------------------------------------------
# percent of cells that express Fbln1

fbln1_exp <- GetAssayData(srat, slot = "data", assay = "RNA") %>% 
  .["Fbln1", ]
length(fbln1_exp[fbln1_exp > 0]) / length(fbln1_exp)

# ------------------------------------------------------------------------------
# make dot plot and violin plots for mouse data

my_VlnPlot <- function(assay, feature) {
  vln <- VlnPlot(srat, assay = assay, features = feature,
                 group.by = "CellType", split.by = "DevStage")
  return(vln)
}

for (assay in c("RNA", "integrated")) {
  DotPlot(srat, assay = assay, cols = c("blue", "blue"),
          features = c("Atxn10", "Fbln1"),
          group.by = "CellType",
          split.by = "DevStage") + 
    scale_y_discrete(limits = rev)
  
  # save dot plot
  ggsave(filename = paste0("expr_dotplot_", assay, ".pdf"), path = date_dir, 
         width = 6, height = 6, units = "in")
  
  my_VlnPlot(assay, "Atxn10") + NoLegend() +
    my_VlnPlot(assay, "Fbln1") +
    plot_layout(ncol = 2)
  
  # save violin plot
  ggsave(filename = paste0("vlnplots_", assay, ".pdf"), path = date_dir,
         width = 10, height = 6, units = "in")
}

# ------------------------------------------------------------------------------
# make new violin plot with P21 before adult

my_VlnPlot("RNA", "Atxn10") + NoLegend() +
  my_VlnPlot("RNA", "Fbln1") +
  plot_layout(ncol = 2)
ggsave(filename = "vlnplots_RNA.pdf", path = date_dir,
       width = 10, height = 6, units = "in")

# ------------------------------------------------------------------------------
# make dot plot for human data (Aldinger et al. 2021)

DotPlot(aldinger_srat, assay = "RNA", cols = rep("blue", 10),
        features = c("ATXN10", "FBLN1"),
        group.by = "figure_clusters",
        split.by = "age") + 
  scale_y_discrete(limits = rev)
ggsave(filename = "human_dotplot_RNA.pdf", path = date_dir,
       width = 8, height = 48, units = "in")

