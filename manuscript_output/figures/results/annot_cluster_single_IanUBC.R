# ==============================================================================
# Annotate tumour cells with SingleR
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(SingleR)

only_80k <- FALSE
srat_rds <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds"
ref_rds <- "/home/rstudio/isilon/private/llau/results/integrated/20240524/25_pc_without_luo/25_pc_rl.qs"
out_dir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/tumour/annot_cluster_single"
UBC_Seurat <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20240825/ubc_subset.qs"

dt <- format(Sys.Date(), "%y%m%d")
out_dir <- sprintf("%s/%s", out_dir, dt)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = FALSE)
}


plotSingleR <- function(tumour_srat, preds, out_dir, cell_order, sfx="",
  predColumn="singleR_preds", pal=NULL, t_pal=NULL) {
  cat("Plotting results...")

  #tumour_srat@meta.data[c("subtype")] <- factor(
  #  tumour_srat@meta.data[c("subtype")],
  #  levels = c("SHH","G3","G4")
  #)

  # prediction heatmap
  cat("heatmap...\n")
  .plt <- plotScoreHeatmap(
    preds,
    annotation_col = tumour_srat@meta.data[c("subtype")],
    cluster_cols = TRUE,
    silent = TRUE
  )
  ggsave(
    sprintf("score_heatmap_%s.pdf", sfx),
    plot = .plt,
    path = out_dir,
    width = 20,
    height = 8,
    units = "in"
  )

cat("barplots...\n")
  # barplots
  .plt1 <- cluster_barplot(
    tumour_srat,
    split.by = predColumn,
    group.by = "subtype",
    position = "fill"
  )

  .plt2 <- cluster_barplot(
    tumour_srat,
    split.by = predColumn,
    group.by = "orig.ident",
    position = "fill"
  ) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  .plt3 <- cluster_barplot(
    tumour_srat,
    split.by = "seurat_clusters",
    group.by = "subtype",
    position = "fill"
  )

  .plt4 <- cluster_barplot(
    tumour_srat,
    split.by = "seurat_clusters",
    group.by = "orig.ident",
    position = "fill"
  ) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  .plt5 <- cluster_barplot(
    tumour_srat,
    split.by = predColumn,
    group.by = "seurat_clusters",
    position = "fill",
    width = 0.5
  )

cat("got past five barplots\n")
  layout <- c(
    area(t = 1, l = 1, r = 1),
    area(t = 1, l = 2, r = 3),
    area(t = 2, l = 1, r = 3)
  )
  .plt <- (
    .plt3 + .plt4 + 
      plot_layout(widths = c(1, 2), guides = "collect") & 
      guides(fill = guide_legend(ncol = 2))
  ) / 
    (.plt1 + .plt2 + .plt5 + plot_layout(design = layout, guides = "collect") & 
      scale_fill_manual(values = pal)
       #scale_fill_manual(values = DiscretePalette(n = length(cell_order), palette = "polychrome"))
     ) + 
    plot_layout(heights = c(1, 2))

  ggsave(
    sprintf("annot_bar_%s.pdf", sfx),
    plot = .plt,
    path = out_dir,
    width = 8,
    height = 8,
    units = "in"
  )

  # rewrite the purrr::map output as a for loop
  .plt <- list()

  for (group.by in c("orig.ident", "subtype", predColumn)) {
    curPal <- t_pal
    if (group.by %in% "orig.ident") {
      curPal <- NULL
    } else if (group.by %in% "subtype") {
      curPal <- t_pal
    } else if (group.by %in% predColumn) {
      curPal <- pal
    }
    .plt[[group.by]] <- DimPlot(
      tumour_srat,
      reduction = "umap",
      group.by = group.by,
      cols = curPal,
      label = FALSE
    )
  }

  ggsave(
    sprintf("umaps_%s.pdf", sfx),
    plot = wrap_plots(.plt),
    path = out_dir,
    width = 16,
    height = 4,
    units = "in"
  )

  cat("done\n")
}

palettes <- list(
   # colour palette for CCA-integrated RL lineage broad cell types
    rl_integ_annot = setNames(
      # remove yellow, too similar to the UBC yellow
      pals::trubetskoy(8)[-3] %>% unname(),
      nm = c("endothelial", "GC", "GCP", "microglia", "oligodendrocyte/OPC", "RL", "UBC")
    ),
    # colour palette for CCA-integrated UBC clusters
    ubc_integ_clust = setNames(
      pals::brewer.set2(6),
      nm = paste0("UBC_", c(0:5))
    ),
    
    mb_subtype = setNames(
      pals::okabe(8)[c(8, 3, 4)],
      nm = c("SHH", "G3", "G4")
    )
)

logFile <- sprintf("%s/log_AldingerSepp_IanUBC.txt", out_dir)
sink(logFile, split=TRUE)
tryCatch({
# ------------------------------------------------------------------------------
# functions

source("/home/rstudio/isilon/private/icheong/CBL_scRNAseq/software/utilities/cluster_barplot.R")

#' Extract normalized gene expression from Seurat object.
#'
#' @param srat Seurat object.
#' @param rerun_norm Rerun `NormalizeData`? Defaults to `TRUE`.
#'
#' @return Matrix of normalized counts.
#'
get_norm_expr <- function(srat, rerun_norm = TRUE) {
  if (rerun_norm) {
    srat <- NormalizeData(srat, assay = "RNA")
  }
  norm_mat <- GetAssayData(srat, slot = "data", assay = "RNA")
  return(norm_mat)
}

# ------------------------------------------------------------------------------
# load Seurat

# tumour 
cat("Reading tumour...")
t0 <- Sys.time()
mb_srat <- readRDS(srat_rds)
print(Sys.time() - t0)
cat("done\n")
mb_srat$subtype <- stringr::str_replace(
  string = mb_srat$orig.ident,
  pattern = "^[:alnum:]*_",
  replacement = ""
)

# reference
cat("Reading reference...")
t0 <- Sys.time()
ref_srat <- qs::qread(ref_rds)
print(Sys.time() - t0)
cat("done\n")

ref_srat <- subset(ref_srat, species == "human")
cat(sprintf("Number of reference cells after subsetting for human: %d\n", ncol(ref_srat)))

ref_srat$common_cell_name <- factor(ref_srat$common_cell_name,
    levels=c("RL","GCP","GC",
        "UBC/GCP progenitor","UBC",
        "microglia",
        "endothelial","oligodendrocyte/OPC"))


# ------------------------------------------------------------------------------
cat("Harmonize cell cluster labels between the UBC set and the RL set")
cat("reading ubc data...")
t0 <- Sys.time()
    UBC_srat <- qs::qread(UBC_Seurat)
    print(Sys.time()-t0)
    cat("done\n")
    
    cat("UBC labels\n")
    ubc_labels <- as.data.frame(as.character(UBC_srat$subclusters))
    ubc_labels$CellID <- rownames(UBC_srat[[]])
    colnames(ubc_labels) <- c("labels","CellID")
    table(ubc_labels$labels,useNA="always")

    cat("RL labels\n")
    RL_labels <- as.data.frame(as.character(ref_srat[[]]$common_cell_name))
    RL_labels$CellID <- rownames(ref_srat[[]])
    colnames(RL_labels) <- c("labels","CellID")
    table(RL_labels$labels,useNA="always")

    # create a new label column that combines the UBC and RL labels using the CellID
    # merge the two dataframes on CellID
    merged_labels <- merge(RL_labels, ubc_labels, by = "CellID", all.x=TRUE)
    idx <- which(!is.na(merged_labels$labels.y))
    #merged_labels$labels.y[idx]<- sprintf("UBC%s",merged_labels$labels.y[idx])
    merged_labels$labels.x[idx] <- merged_labels$labels.y[idx]

    midx <- match(rownames(ref_srat[[]]),merged_labels$CellID)
    if (all.equal(merged_labels$CellID[midx],rownames(ref_srat[[]]))!= TRUE) {
        stop("CellID mismatch")
    }

ref_srat[[]]$RLwithUBC <- merged_labels$labels.x[midx]

# remove cells with RLwithUBC label UBC
cat("Removing cells with RLwithUBC label UBC\n")
keep_idx <- which(!ref_srat[[]]$RLwithUBC %in% c("UBC"))
ref_srat <- subset(ref_srat, cells = Cells(ref_srat)[keep_idx])    

# relabel all "GN" to "GC" for consistency
cat("Relabelling all GN to GC\n")
gn_idx <- which(ref_srat$common_cell_name == "GN")
ref_srat$common_cell_name[gn_idx] <- "GC"

ref_srat$RLwithUBC <- factor(ref_srat$RLwithUBC,
    levels=c("RL","GCP","GC",
        "UBC", "UBC_0","UBC_1","UBC_2","UBC_3","UBC_4","UBC_5",
        "microglia",
        "endothelial","oligodendrocyte/OPC"))


#ref_srat$common_cell_name <- ref_srat$RLwithUBC

# ------------------------------------------------------------------------------
# SingleR predictions

# re-run normalization (SingleR needs in normalized counts, don't use
# SCTransform) and subset the normalized expression matrices
cat("Normalizing tumour...")
mb_mat <- get_norm_expr(mb_srat, rerun_norm = TRUE)
cat("done\n")
cat("Normalizing reference...")
ref_mat <- get_norm_expr(ref_srat, rerun_norm = TRUE)
cat("done\n")



# note de.method is "wilcox" cuz the reference is a single cell dataset
# (see https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html)
message("Running SingleR")
t0 <- Sys.time()
annot_prediction <- SingleR(
  test = mb_mat,
  ref = ref_mat,
  labels = ref_srat$common_cell_name,
  de.method = "wilcox",
  num.threads = 24
)
print(Sys.time() - t0)

# factor the label
cell_order <- levels(ref_srat$common_cell_name)
annot_prediction$labels <- factor(annot_prediction$labels, levels = cell_order)

# add labels back to metadata
mb_srat$singleR_preds <- annot_prediction[rownames(mb_srat[[]]), "labels"]
plotSingleR(mb_srat, annot_prediction, 
  sfx="AldingerSepp",
  out_dir=out_dir,
  cell_order=cell_order,
  predColumn="singleR_preds",
  pal = palettes$rl_integ_annot,
  t_pal = palettes$mb_subtype
  )

orig_mb_srat <- mb_srat
orig_ref_srat <- ref_srat

# calculate and print out what proportion of each subtype is assigned to each cell type
cat("Calculating proportions of cell types in each MB subtype...\n")
prop_table <- table(orig_mb_srat$subtype, orig_mb_srat$singleR_preds)
prop_table <- (prop_table / rowSums(prop_table))*100
print(round(prop_table,2))

# ------------------------------------------------------------------------------
cat("now we subset for UBC-like cells only in Group 3 and 4 tumours and rerun SingleR but with the RLwithUBC labels\n")
# subset to only Group 3 and 4 tumours and UBC-like cells
ubc_like_idx <- which(
  mb_srat$subtype %in% c("G3","G4") &
    mb_srat$singleR_preds %in% "UBC"
)

mb_srat <- subset(mb_srat, cells = Cells(mb_srat)[ubc_like_idx])
cat(sprintf("Number of UBC-like cells in Group 3 and 4 tumours: %d\n", ncol(mb_srat)))

# now subset the reference for just UBC_0 to UBC_5
ubc_ref_idx <- which(
  ref_srat$RLwithUBC %in% c("UBC_0","UBC_1","UBC_2","UBC_3","UBC_4","UBC_5")
)
ref_srat <- subset(ref_srat, cells = Cells(ref_srat)[ubc_ref_idx])
cat(sprintf("Number of UBC reference cells: %d\n", ncol(ref_srat)))

# re-run normalization (SingleR needs in normalized counts, don't use
# SCTransform) and subset the normalized expression matrices
cat("Normalizing tumour...")
mb_mat <- get_norm_expr(mb_srat, rerun_norm = TRUE)
cat("done\n")
cat("Normalizing reference...")
ref_mat <- get_norm_expr(ref_srat, rerun_norm = TRUE)
cat("done\n")


# note de.method is "wilcox" cuz the reference is a single cell dataset
# (see https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html)
message("Running SingleR")
t0 <- Sys.time()
annot_prediction_UBC <- SingleR(
  test = mb_mat,
  ref = ref_mat,
  labels = ref_srat$RLwithUBC,
  de.method = "wilcox",
  num.threads = 24
)
print(Sys.time() - t0)  

# factor the label
cell_order <- paste("UBC",0:5,sep="_")
annot_prediction_UBC$labels <- factor(annot_prediction_UBC$labels, 
  levels = cell_order)
# add labels back to metadata
mb_srat$singleR_UBCsubtype <- annot_prediction_UBC[rownames(mb_srat[[]]), "labels"]
plotSingleR(mb_srat, annot_prediction_UBC, sfx="AldingerSepp_UBCsubtypes", 
  out_dir=out_dir,
  cell_order=cell_order,
  predColumn="singleR_UBCsubtype",
  pal = palettes$ubc_integ_clust,
  t_pal = palettes$mb_subtype
)

cat("Calculating proportions of UBC subtypes in each MB subtype...\n")
prop_table <- table(mb_srat$subtype, mb_srat$singleR_UBCsubtype)
prop_table <- (prop_table / rowSums(prop_table))*100
print(round(prop_table,2))


}, error=function(ex){
    print(ex)
}, finally={
    sink()
    message("\n***SESSION INFO***\n")
print(sessionInfo())
})


# ------------------------------------------------------------------------------



