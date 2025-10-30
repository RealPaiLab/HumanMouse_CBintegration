# ==============================================================================
# Integrate RL-derived Vladoiu mouse cells with Liam's RL-derived human cells
# using harmony.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/CBL_scRNAseq/software/")

library(biomaRt)
library(tidyverse)
library(Seurat)
library(harmony)

# set data and output directories
root_dir <- "CBL_scRNAseq"
human_data_dir <- file.path("/isilon", root_dir, "data/human")
mouse_data_dir <- file.path("/isilon", root_dir, "data/mouse")

human_out_dir <- file.path("", root_dir, "results/human")
mouse_out_dir <- file.path("", root_dir, "results/mouse")
integ_out_dir <- file.path(paste0("/", root_dir), "results/integrated")

# integ_date_dir <- file.path(integ_out_dir, format(Sys.Date(), "%Y%m%d"))
integ_date_dir <- file.path(integ_out_dir, "20221004")

if (!dir.exists(integ_date_dir)) {
  dir.create(integ_date_dir)
}

# ------------------------------------------------------------------------------
# import data

liam_srat <- readRDS(file.path(human_data_dir, "Aldinger/glutamatergic_dev_Liam.RDS"))
vladRL_srat <- readRDS(file.path(mouse_out_dir, "Vladoiu/merged_seurat_RLonly.rds"))

# ------------------------------------------------------------------------------
# convert mouse genes to human genes

# if file of orthologous genes exists, use that one; otherwise fetch from biomart
if (file.exists(file.path(integ_out_dir, "hgnc_mgi_orth_genes.csv"))) {
  m2h_genes <- read.csv(file.path(integ_out_dir, "hgnc_mgi_orth_genes.csv"))
} else {
  # load biomart/ensembl databases
  source("./utilities/convert_genes.R")
  human_bm <- load_human_bm()
  mouse_bm <- load_mouse_bm()
  
  # get mouse genes and convert to human genes
  mouse_genes <- rownames(vladRL_srat)
  m2h_genes <- convert_m2h(mouse_genes, mouse_bm, human_bm)
  
  # save the orthologous genes to a CSV file
  # readr::write_csv(x = m2h_genes, file = file.path(integ_out_dir, "vladoiu_orth_genes.csv"))
}

# remove all duplicated genes
# also remove Pisd (some wonky stuff is happening with that gene)
m2h_genes <- m2h_genes[duplicated(m2h_genes$MGI.symbol) == FALSE
                       & duplicated(m2h_genes$HGNC.symbol) == FALSE
                       & m2h_genes$MGI.symbol != "Pisd", ]

# rename mouse genes to their orthologous human genes
source(file.path("", root_dir, "software/utilities/convert_genes.R"))
vlad_srat_hgene <- rename_genes(vladRL_srat, old_genes = m2h_genes$MGI.symbol, new_genes = m2h_genes$HGNC.symbol)

# ------------------------------------------------------------------------------
# integrate datasets using harmony

# first step is to merge the human and mouse into a single Seurat object
merge_srat <- merge(liam_srat, vlad_srat_hgene)

# add species to metadata
merge_srat@meta.data <- merge_srat@meta.data %>% 
  mutate(
    dataset = as.factor(case_when(
      orig.ident == "80k" ~ "Aldinger",
      orig.ident != "80k" ~ "Vladoiu"
    )),
    species = as.factor(case_when(
      orig.ident == "80k" ~ "human",
      orig.ident != "80k" ~ "mouse"
    ))
  )

# set variable features to common genes
VariableFeatures(merge_srat) <- SelectIntegrationFeatures(
  list(liam_srat, vlad_srat_hgene),
  nfeatures = 5000
)

# run PCA
merge_srat <- RunPCA(merge_srat, features = VariableFeatures(merge_srat), npcs = 100)

# run harmony
integ_srat <- RunHarmony(
  object = merge_srat,
  group.by.vars = "species",
  plot_convergence = TRUE,
  assay.use = "SCT",
  project.dim = FALSE
)

# ------------------------------------------------------------------------------
# clustering and dimensional reduction

# integ_srat <- RunPCA(integ_srat, features = VariableFeatures(integ_srat), npcs = 100)
# DimHeatmap(integ_srat, dims = 1:15, cells = 500)
ElbowPlot(integ_srat, ndims = 50, reduction = "harmony")

# set ndims for clustering and dimensional reductions (based on elbow plot)
ndims <- 25

# cluster cells
integ_srat <- FindNeighbors(integ_srat, reduction = "harmony", dims = 1:ndims) %>% 
  FindClusters(.)

# run t-SNE
integ_srat <- RunTSNE(integ_srat, reduction = "harmony", dims = 1:ndims, num_threads = 32)

# run UMAP
integ_srat <- RunUMAP(integ_srat, reduction = "harmony", dims = 1:ndims)

# ------------------------------------------------------------------------------
# make and save plots

# make tsne and umap plots, grouped by seurat clusters and by species
for (reduc in c("tsne", "umap")) {
  
  # seurat clusters
  plt1 <- DimPlot(
    integ_srat,
    reduction = reduc,
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  ) + 
    NoLegend()
  
  # species
  plt2 <- DimPlot(
    integ_srat,
    reduction = reduc,
    group.by = "species",
    label = FALSE
  )
  
  # combine plots and save
  fname <- paste0(reduc, ".png")
  ggsave(
    filename = fname,
    plot = plt1 + plt2,
    path = integ_date_dir,
    width = 16,
    height = 8,
    units = "in",
    dpi = 600
  )
}

# make umap plots grouped by cell types
source(file.path("", root_dir, "software/mouse/Vladoiu/add_annotations.R"))
integ_srat <- label_vladoiu_cells(integ_srat)

umap_human_cells <- DimPlot(
  integ_srat,
  cells = rownames(integ_srat@meta.data[integ_srat$species == "human", ]),
  reduction = "umap",
  group.by = "new_cell_type",
  label = TRUE,
  repel = TRUE
) + 
  ggtitle("human cell types") + 
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2)))

umap_mouse_cells <- DimPlot(
  integ_srat,
  cells = rownames(integ_srat@meta.data[integ_srat$species == "mouse", ]),
  cols = c(
    # RL/UBC lineage
    scales::brewer_pal(palette = "Reds")(4),
    # granule cell (progenitors)
    scales::brewer_pal(palette = "Blues")(3),
    # other glutamatergic
    scales::brewer_pal(palette = "YlOrBr")(3),
    # stem cells
    scales::brewer_pal(palette = "BuPu")(2),
    # GABAergic lineage
    scales::brewer_pal(palette = "Greens")(5),
    # glial/non-neuronal cells
    scales::brewer_pal(palette = "PuRd")(5),
    # NA
    "#DFDFDF"
  ),
  reduction = "umap",
  group.by = "mouse_cell_type",
  label = FALSE
) + 
  ggtitle("mouse cell types") + 
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))

ggsave(
  filename = "umap_celltypes.png",
  plot = umap_human_cells + umap_mouse_cells,
  path = integ_date_dir,
  width = 15,
  height = 10,
  units = "in",
  dpi = 600
)

saveRDS(integ_srat, file = file.path(integ_date_dir, "vladoiu_liam_RL_harmony.rds"))
