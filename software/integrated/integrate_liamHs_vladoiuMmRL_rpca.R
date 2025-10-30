# ==============================================================================
# Integrate RL-derived Vladoiu mouse cells with Liam's RL-derived human cells
# using Seurat's RPCA method.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/CBL_scRNAseq/software/")

library(biomaRt)
library(tidyverse)
library(Seurat)

# set data and output directories
root_dir <- "CBL_scRNAseq"
human_data_dir <- file.path("/isilon", root_dir, "data/human")
mouse_data_dir <- file.path("/isilon", root_dir, "data/mouse")

human_out_dir <- file.path("", root_dir, "results/human")
mouse_out_dir <- file.path("", root_dir, "results/mouse")
integ_out_dir <- file.path(paste0("/", root_dir), "results/integrated")

integ_date_dir <- file.path(integ_out_dir, format(Sys.Date(), "%Y%m%d"))

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
# integrate datasets using RPCA
# (https://satijalab.org/seurat/articles/integration_rpca.html)

# Make sure that SCTransform has already been called on the Seurat objects, i.e.
# object@active.assay should be "SCT". This code assumes that the Seurat objects
# have been separately normalized with SCTransform.

# combine human and mouse Seurat into list
obj_list <- list(liam_srat, vlad_srat_hgene)

# prep for integration
obj_list <- lapply(
  X = obj_list,
  FUN = SCTransform,
  variable.features.n = 5000,
  vars.to.regress = "CC.Difference",
  return.only.var.genes = FALSE
)
features <- SelectIntegrationFeatures(obj_list, nfeatures = 5000)
obj_list <- lapply(
  X = obj_list,
  FUN = RunPCA,
  features = features
)
obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)

# iterate through different number of k.anchor for FindIntegrationAnchors
k.anchors <- 5*2^seq(0, 6)
for (k in k.anchors) {
  message(sprintf("Integrating with k.anchors = %s", k))
  
  # find anchors and integrate
  anchors <- FindIntegrationAnchors(
    obj_list,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = "rpca",
    k.anchor = k
  )
  integ_srat <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  
  # run PCA, check PCs
  integ_srat <- RunPCA(integ_srat, features = VariableFeatures(integ_srat), npcs = 100)
  DimHeatmap(integ_srat, dims = 1:15, cells = 500)
  ElbowPlot(integ_srat, ndims = 100)
  
  # ---
  # clustering and dimensional reduction
  
  # set ndims for clustering and dimensional reductions (based on elbow plot)
  ndims <- 30
  
  # cluster cells
  integ_srat <- FindNeighbors(integ_srat, dims = 1:ndims) %>% 
    FindClusters(.)
  
  # run t-SNE
  integ_srat <- RunTSNE(integ_srat, dims = 1:ndims, num_threads = 32)
  
  # run UMAP
  integ_srat <- RunUMAP(integ_srat, dims = 1:ndims)
  
  # ---
  # add species metadata
  
  integ_srat@meta.data <- integ_srat@meta.data %>% 
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
  
  # ---
  # save plots and Seurat object
  
  # create directory
  k_path <- file.path(integ_date_dir, paste0("kanchors", k))
  if (!dir.exists(k_path)) {
    dir.create(k_path)
  }
  
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
      path = k_path,
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
    guides(colour = guide_legend(nrow = 1))
  
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
    path = k_path,
    width = 15,
    height = 10,
    units = "in",
    dpi = 600
  )
  
  saveRDS(integ_srat, file = file.path(k_path, "vladoiu_liam_RL_rpca.rds"))
}
