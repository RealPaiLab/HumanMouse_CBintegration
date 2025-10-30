# ==============================================================================
# Annotate tumour cells with SingleR
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(SingleR)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # Seurat tumour file
  "--srat_rds",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # Seurat object for reference
  "--ref_rds",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # only use the 80k Aldinger dataset
  "--only_80k",
  action = "store_true"
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  # for testing and troubleshooting
  args <- parser$parse_args(c(
    "--srat_rds", "/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds",
    "--ref_rds", "/isilon/CBL_scRNAseq-archived/data/human/Aldinger/seurat.rds",
    "--out_dir", "/CBL_scRNAseq/results/tumour/Vladoiu/20230601/"
  ))
} else {
  args <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(args$out_dir)) {
  dir.create(args$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# functions

source("/CBL_scRNAseq/software/utilities/cluster_barplot.R")

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
mb_srat <- readRDS(args$srat_rds)
mb_srat$subtype <- stringr::str_replace(
  string = mb_srat$orig.ident,
  pattern = "^[:alnum:]*_",
  replacement = ""
)

# reference
ref_srat <- readRDS(args$ref_rds)

# only use Aldinger 80k for reference
if (args$only_80k) {
  ref_srat <- subset(ref_srat, subset = orig.ident == "80k")
}

# ------------------------------------------------------------------------------
# copy Hendrikse annotations to Aldinger metadata

rl_celltypes <- readRDS("/isilon/CBL_scRNAseq-archived/data/human/Aldinger/glutamatergic_dev_Liam.RDS") %>% 
  .@meta.data %>% 
  select(new_cell_type)

ref_srat@meta.data <- merge(
  x = ref_srat@meta.data,
  y = rl_celltypes,
  by = 0, # 0 means the row names
  all.x = TRUE,
  sort = FALSE
) %>% 
  mutate(
    merged_cell_type = case_when(
      !is.na(new_cell_type) ~ as.character(new_cell_type),
      TRUE ~ as.character(figure_clusters)
    )
  ) %>% 
  column_to_rownames("Row.names")

# ------------------------------------------------------------------------------
# SingleR predictions

# re-run normalization (SingleR needs in normalized counts, don't use
# SCTransform) and subset the normalized expression matrices
mb_mat <- get_norm_expr(mb_srat, rerun_norm = TRUE)
ref_mat <- get_norm_expr(ref_srat, rerun_norm = TRUE)

# get reference annotations (make sure it's in the same order as expr matrix)
ref_labels <- ref_srat@meta.data[colnames(ref_mat), "merged_cell_type"]

# note de.method is "wilcox" cuz the reference is a single cell dataset
# (see https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html)
message("Running SingleR")
annot_prediction <- SingleR(
  test = mb_mat,
  ref = ref_mat,
  labels = ref_labels,
  de.method = "wilcox"
)

# factor the label
cell_order <- c(
  "RL-VZ",
  "RL-SVZ",
  "Early UBCs",
  "Late UBCs",
  "GCP",
  "Early GN",
  "GN",
  levels(ref_srat$figure_clusters)
)
annot_prediction$labels <- factor(annot_prediction$labels, levels = cell_order)

# add labels back to metadata
mb_srat$singleR_preds <- annot_prediction[rownames(mb_srat[[]]), "labels"]

# ------------------------------------------------------------------------------
# visualizations

# prediction heatmap
.plt <- plotScoreHeatmap(
  annot_prediction,
  annotation_col = mb_srat@meta.data[c("subtype")],
  cluster_cols = TRUE,
  silent = TRUE
)
ggsave(
  "score_heatmap.png",
  plot = .plt,
  path = args$out_dir,
  width = 20,
  height = 8,
  units = "in",
  dpi = 600
)


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
  path = args$out_dir,
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
  path = args$out_dir,
  width = 16,
  height = 4,
  units = "in",
  dpi = 600
)


# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

