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
ref_rds <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/data/human/Aldinger/glutamatergic_dev_Liam.RDS"
out_dir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/tumour/annot_cluster_single"

dt <- format(Sys.Date(), "%y%m%d")
out_dir <- sprintf("%s/%s", out_dir, dt)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = FALSE)
}

logFile <- sprintf("%s/log.txt", out_dir)
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
ref_srat <- readRDS(ref_rds)
print(Sys.time() - t0)
cat("done\n")

# only use Aldinger 80k for reference
if (only_80k) {
  ref_srat <- subset(ref_srat, subset = orig.ident == "80k")
}

# ------------------------------------------------------------------------------
# copy Hendrikse annotations to Aldinger metadata
###cat("Merging cell type annotations...")
###rl_celltypes <- readRDS(liamFile) %>% 
###  .@meta.data %>% 
###  select(new_cell_type)
###
###ref_srat@meta.data <- merge(
###  x = ref_srat@meta.data,
###  y = rl_celltypes,
###  by = 0, # 0 means the row names
###  all.x = TRUE,
###  sort = FALSE
###) %>% 
###  mutate(
###    merged_cell_type = case_when(
###      !is.na(new_cell_type) ~ as.character(new_cell_type),
###      TRUE ~ as.character(figure_clusters)
###    )
###  ) %>% 
###  column_to_rownames("Row.names")
###cat("done\n")

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
annot_prediction <- SingleR(
  test = mb_mat,
  ref = ref_mat,
  labels = ref_srat$new_cell_type,
  de.method = "wilcox"
)

# factor the label
cell_order <- levels(ref_srat$new_cell_type)
annot_prediction$labels <- factor(annot_prediction$labels, levels = cell_order)

# add labels back to metadata
mb_srat$singleR_preds <- annot_prediction[rownames(mb_srat[[]]), "labels"]

# ------------------------------------------------------------------------------
# visualizations

cat("Plotting results...")
###
#### prediction heatmap
###.plt <- plotScoreHeatmap(
###  annot_prediction,
###  annotation_col = mb_srat@meta.data[c("subtype")],
###  cluster_cols = TRUE,
###  silent = TRUE
###)
###ggsave(
###  "score_heatmap.png",
###  plot = .plt,
###  path = out_dir,
###  width = 20,
###  height = 8,
###  units = "in",
###  dpi = 600
###)


# barplots
.plt1 <- cluster_barplot(
  mb_srat,
  split.by = "singleR_preds",
  group.by = "subtype",
  position = "fill"
)

.plt2 <- cluster_barplot(
  mb_srat,
  split.by = "singleR_preds",
  group.by = "orig.ident",
  position = "fill"
) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

.plt3 <- cluster_barplot(
  mb_srat,
  split.by = "seurat_clusters",
  group.by = "subtype",
  position = "fill"
)

.plt4 <- cluster_barplot(
  mb_srat,
  split.by = "seurat_clusters",
  group.by = "orig.ident",
  position = "fill"
) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

.plt5 <- cluster_barplot(
  mb_srat,
  split.by = "singleR_preds",
  group.by = "seurat_clusters",
  position = "fill",
  width = 0.5
)

layout <- c(
  area(t = 1, l = 1, r = 1),
  area(t = 1, l = 2, r = 3),
  area(t = 2, l = 1, r = 3)
)
.plt <- (
  .plt3 + .plt4 + 
    plot_layout(widths = c(1, 2), guides = "collect") & 
    guides(fill = guide_legend(ncol = 2))
    # scale_fill_manual(values = DiscretePalette(
    #   length(unique(mb_srat$seurat_clusters)),
    #   "stepped"
    # ))
) / 
  (.plt1 + .plt2 + .plt5 + plot_layout(design = layout, guides = "collect") & 
     scale_fill_manual(values = DiscretePalette(n = length(cell_order), palette = "polychrome"))
   ) + 
  plot_layout(heights = c(1, 2))

ggsave(
  "annot_bar.png",
  plot = .plt,
  path = out_dir,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)


# UMAPs
.plt <- purrr::map(
  .x = c("orig.ident", "subtype", "singleR_preds"),
  .f = \(group.by) {
    .plt <- DimPlot(
      mb_srat,
      reduction = "umap",
      group.by = group.by,
      label = FALSE
    )
    return(.plt)
  }
)

ggsave(
  "umaps.png",
  plot = wrap_plots(.plt),
  path = out_dir,
  width = 16,
  height = 4,
  units = "in",
  dpi = 600
)

cat("done\n")

}, error=function(ex){
    print(ex)
}, finally={
    sink()
 #   message("\n***SESSION INFO***\n")
#print(sessionInfo())
})


# ------------------------------------------------------------------------------



