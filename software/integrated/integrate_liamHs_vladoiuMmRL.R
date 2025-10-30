# ==============================================================================
# Integrate RL-derived Vladoiu mouse cells with Liam's RL-derived human cells
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

# ==============================================================================
# import data
# ==============================================================================

liam_srat <- readRDS(file.path(human_data_dir, "Aldinger/glutamatergic_dev_Liam.RDS"))
vladRL_srat <- readRDS(file.path(mouse_out_dir, "Vladoiu/merged_seurat_RLonly.rds"))

# ==============================================================================
# convert mouse genes to human genes
# ==============================================================================

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

# change mouse gene names to human gene names
# if human ortholog does not exist, keep mouse gene name
get_new_names <- function(
  old_genes,
  new_genes,
  old_names
) {
  # convert old gene names to new gene names using the genes passed into the function
  new_names <- sapply(
    X = old_names, 
    FUN = function(X) {
      if (X %in% old_genes) {
        # if an orthologous gene exists, then set the name to the new gene
        X <- new_genes[old_genes == X]
      } else {
        # if no orthologous gene was found, keep the old (original) gene name
        X
      }
    }, 
    USE.NAMES = FALSE
  )
  
  return(new_names)
}
rename_genes <- function(
  object,
  old_genes,
  new_genes
) {
  # loop through assays
  for (assay in names(object@assays)) {
    
    # loop through counts, data, and scale.data slots
    for (slot in c("counts", "data", "scale.data")) {
      
      # get vector of old gene names
      old_names <- rownames(slot(object[[assay]], slot))
      
      # continue to next loop iteration if slot is empty
      if (length(old_names) == 0) {next}
      
      # get the new gene names
      new_names <- get_new_names(old_genes, new_genes, old_names)
      
      # set the new gene names
      rownames(slot(object[[assay]], slot)) <- unlist(new_names)
    }
    
    # set new names for var.features slot
    if (length(object[[assay]]@var.features) == 0) {
      next
    }
    old_names <- object[[assay]]@var.features
    object[[assay]]@var.features <- unlist(get_new_names(old_genes, new_genes, old_names))
  }
  
  return(object)
}

vlad_srat_hgene <- rename_genes(vladRL_srat, old_genes = m2h_genes$MGI.symbol, new_genes = m2h_genes$HGNC.symbol)

# ==============================================================================
# integrate the datasets
# ==============================================================================

obj_list <- list(liam_srat, vlad_srat_hgene)

# re-run sctransform
obj_list <- lapply(X = obj_list, FUN = SCTransform,
                   variable.features.n = 5000,
                   vars.to.regress = "CC.Difference",
                   return.only.var.genes = FALSE)

# prep for integration
features <- SelectIntegrationFeatures(obj_list, nfeatures = 5000)
obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)

# find anchors and integrate
anchors <- FindIntegrationAnchors(obj_list, normalization.method = "SCT", 
                                  anchor.features = features, reduction = "cca")
integ_srat <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# ==============================================================================
# clustering and dimensional reduction
# ==============================================================================

# run PCA, check PCs
integ_srat <- RunPCA(integ_srat, features = VariableFeatures(integ_srat), npcs = 100)
DimHeatmap(integ_srat, dims = 1:15, cells = 500)
ElbowPlot(integ_srat, ndims = 100)

# set ndims for clustering and dimensional reductions (based on elbow plot)
ndims <- 30

# cluster cells
integ_srat <- FindNeighbors(integ_srat, dims = 1:ndims) %>% 
  FindClusters(.)

# run t-SNE
integ_srat <- RunTSNE(integ_srat, dims = 1:ndims, num_threads = 32)

# run UMAP
integ_srat <- RunUMAP(integ_srat, dims = 1:ndims)

# ==============================================================================
# save integrated Seurat plots and object
# ==============================================================================

# add dataset and species to metadata
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

# save plots
files <- c("integ_tsne_clusters.pdf",
           "integ_tsne_species.pdf",
           "integ_umap_clusters.pdf",
           "integ_umap_species.pdf")
i <- 1
for (reduc in c("tsne", "umap")) {
  for (group in c("seurat_clusters", "species")) {
    plt <- DimPlot(integ_srat, reduction = reduc, label = TRUE, repel = TRUE, group.by = group)
    ggsave(filename = files[i], plot = plt, path = integ_date_dir, 
           width = 10, height = 7.5, units = "in")
    i <- i+1
  }
}

# save object
saveRDS(object = integ_srat,
        file = file.path(integ_out_dir, "vladoiu_liam_RL.rds"))
