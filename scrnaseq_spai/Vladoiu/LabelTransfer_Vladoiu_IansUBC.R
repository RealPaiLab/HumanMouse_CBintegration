rm(list=ls())

require(Seurat)
require(ggplot2)
require(ggpubr)
require(SingleR)
require(SingleCellExperiment)
require(scran)
library(dplyr)

mb <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds"

devFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/assignClusterIdentity/cca_RLlineage_only_241107.qs"

#UBC_Seurat <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"
UBC_Seurat <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20240825/ubc_subset.qs"
#UBCclusters <- "SCT_snn_res.0.5"

  cat("Loading Vladoiu data\n")
    t0 <- Sys.time()
    tumour_srat <- readRDS(mb)
    print(Sys.time()-t0)


 print(sprintf("Tumour data: %i samples, %i cells, %i features",
        length(unique(tumour_srat$orig.ident)),
        ncol(tumour_srat),
        nrow(tumour_srat)))
    cat("\n")
    
   cat("SingleR projections not found, running SingleR\n")
    cat("Now loading Aldinger/Sepp dev data\n")
    t0 <- Sys.time()
    rl_srat <- readRDS(devFile)
    print(Sys.time()-t0)
    # count number of cells and genes in rl_srat
    cat("Number of cells and genes in rl_srat\n")
    cat(sprintf("Cells: %i\n", ncol(rl_srat)))
    cat(sprintf("Genes: %i\n", nrow(rl_srat)))
    
    cat("Now loading UBC clusters\n")
    t0 <- Sys.time()
    cat("\n")

    cat("Now loading UBC clusters\n")
    #UBC_srat <- readRDS(UBC_Seurat)
    UBC_srat <- qs::qread(UBC_Seurat)
    print(Sys.time()-t0)
    cat("Number of cells and genes in UBC_srat\n")
    cat(sprintf("Cells: %i\n", ncol(UBC_srat)))
    cat(sprintf("Genes: %i\n", nrow(UBC_srat)))
    cat("Normalizing data\n")
    ####UBC_mat <- get_norm_expr(UBC_srat, rerun_norm = TRUE)
    cat("\n")
    
    DefaultAssay(rl_srat) <- "integrated"    
    DefaultAssay(UBC_srat) <- "integrated"

    anchors <- FindTransferAnchors(
      reference = rl_srat,
      query = tumour_srat)

predictions <- TransferData(anchorset = anchors,
    refdata = rl_srat$snn_res.0.8_celltype)

tumour_srat <- AddMetaData(tumour_srat, metadata = predictions)

# print predictions table for each tumour subtype
print(table(tumour_srat$predicted.id, tumour_srat$subtype))