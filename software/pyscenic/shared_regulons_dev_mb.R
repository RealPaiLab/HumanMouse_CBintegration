# ==============================================================================
# Find common regulons for cell types between human cerebellar development and
# medulloblastoma.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--out_dir", "/CBL_scRNAseq/results/pyscenic/20240223"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# load functions
source("/CBL_scRNAseq/software/utilities/pyscenic_regulons.R")

# ------------------------------------------------------------------------------
# load data

# regulons and expression
dev_reg <- read.csv("/CBL_scRNAseq/results/pyscenic/20230725/integ/aldinger_RL.rss.csv", row.names = 1)
mb_reg <- read.csv("/CBL_scRNAseq/results/pyscenic/20231018/vladoiu_mb.rss.csv", row.names = 1)

# MB subtypes
mb_srat <- read_rds("/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds")

# ------------------------------------------------------------------------------
# get top n regulons for each cell type

top_n_regs <- 8

# >>> developing cerebellum

# take cell types with large proportion of human cells (see 2022-09-11 bar plots)
dev_cell_types <- c(
  "5-Human GN and Mouse UBC",
  "6-Human GN",
  "7-Homol UBC",
  "8-Human GCP/RL_SVZ",
  "10-Human GN/Mouse UBC",
  "12-Human RL_VZ and Mouse GCP/progenitor",
  "19-NonHomol UBC",
  "20-NonHomol UBC",
  "21-Human GN"
)

top_dev_regs <- map(
  .x = dev_cell_types,
  .f = \(x, rss = dev_reg, num_regulons = top_n_regs) {
    # function from `pyscenic_regulons.R`
    top_regs <- get_top_regulons(x, rss = rss, num_regulons = num_regulons)
    return(top_regs)
  }
) %>% 
  `names<-`(dev_cell_types)

write_rds(
  x = top_dev_regs,
  file = file.path(arg_list$out_dir, "top_dev_regulons.rds")
)

# <<<

# >>> medulloblastoma

# sort MB cell clusters
mb_cell_types <- rownames(mb_reg) %>% str_sort(numeric = TRUE)

top_mb_regs <- map(
  .x = mb_cell_types,
  .f = \(x, rss = mb_reg, num_regulons = top_n_regs) {
    # function from `pyscenic_regulons.R`
    top_regs <- get_top_regulons(x, rss = rss, num_regulons = num_regulons)
    return(top_regs)
  }
) %>% 
  `names<-`(mb_cell_types)

write_rds(
  x = top_mb_regs,
  file = file.path(arg_list$out_dir, "top_mb_regulons.rds")
)

# <<<

# ------------------------------------------------------------------------------
# make matrix for heatmap

num_shared_regs <- expand_grid(dev_cell_types, mb_cell_types) %>% 
  rowwise() %>% # rowwise is necessary
  mutate(
    # get list of shared shared regulons
    intersection = list(intersect(
      top_dev_regs[[dev_cell_types]],
      top_mb_regs[[mb_cell_types]]
    )),
    # get number of shared regulons
    intersection_size = length(intersection),
    # collapse list of shared regulons into single string
    intersection = case_when(
      intersection_size > 0 ~ paste(intersection, collapse = ", "),
      intersection_size == 0 ~ NA_character_
    )
  )

message("Saving shared regulons")
write_csv(x = num_shared_regs, file = file.path(arg_list$out_dir, "shared_regulons.csv"))

mat <- num_shared_regs %>% 
  select(-intersection) %>% 
  pivot_wider(
    names_from = mb_cell_types,
    values_from = intersection_size
  ) %>% 
  column_to_rownames("dev_cell_types") %>% 
  as.matrix()

message("Saving matrix for heatmap")
write.csv(
  x = mat,
  file = file.path(arg_list$out_dir, "shared_regulons_mat.csv"),
  row.names = TRUE
)
saveRDS(mat, file = file.path(arg_list$out_dir, "shared_regulons_mat.rds"))

# ------------------------------------------------------------------------------
# make subtype bar plot annotation

subtype_mat <- table(mb_srat$seurat_clusters, mb_srat$subtype)

# reorder so SHH comes first
subtype_mat <- subtype_mat[, c("SHH", "G3", "G4")]

# scale per cluster (so there is a proportion of subtypes for each cluster)
subtype_mat <- scale(
  t(subtype_mat),
  center = FALSE,
  scale = rowSums(subtype_mat)
) %>% 
  t()

# bar plot colours
subtype_cols <- pals::brewer.set2(5)[3:5]

subtype_annot <- HeatmapAnnotation(
  subtype = anno_barplot(subtype_mat,
                         height = unit(0.5, "in"),
                         gp = gpar(fill = subtype_cols)),
  show_legend = TRUE,
  annotation_label = c(subtype = "MB\nsubtype")
)

# make legend
lgd_list <- list(
  Legend(
    labels = colnames(subtype_mat),
    title = "MB\nsubtype",
    legend_gp = gpar(fill = subtype_cols)
  )
)

# save matrix
message("Saving matrix for subtype bar plot")
write.csv(
  x = subtype_mat,
  file = file.path(arg_list$out_dir, "subtype_annot.csv"),
  row.names = TRUE
)
saveRDS(subtype_mat, file = file.path(arg_list$out_dir, "subtype_annot.rds"))

# ------------------------------------------------------------------------------
# make heatmap

# colour function for heatmap
col_fun <- circlize::colorRamp2(
  breaks = c(0, max(mat)),
  hcl_palette = "Viridis"
)

# row/column title format
title_gp <- gpar(fontface = "bold")

# row/column names format
names_gp <- gpar(fontsize = 10)

# row labels
row_labels <- structure(
  rownames(mat),
  names = rownames(mat)
)
row_labels["7-Homol UBC"] <- "7-Common UBC"
row_labels["19-NonHomol UBC"] <- "19-Human-specific UBC"
row_labels["20-NonHomol UBC"] <- "20-Human-specific UBC"

hm <- Heatmap(
  matrix = mat,
  col = col_fun,
  name = "shared\nregulons",
  rect_gp = gpar(col = "white", lwd = 2),
  row_title = "developing cerebellum\ncell types",
  row_title_side = "left",
  row_title_gp = title_gp,
  column_title = "MB cell clusters",
  column_title_side = "bottom",
  column_title_gp = title_gp,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_labels = row_labels[rownames(mat)],
  row_names_side = "left",
  row_names_max_width = max_text_width(rownames(mat)),
  row_names_gp = names_gp,
  column_names_rot = 0,
  column_names_gp = names_gp,
  column_names_centered = TRUE,
  top_annotation = subtype_annot,
  width = ncol(mat) * unit(0.25, "in"),
  height = nrow(mat) * unit(0.25, "in")
)
# need `draw()` to calculate heatmap width/height
hm <- draw(hm, annotation_legend_list = lgd_list, merge_legend = TRUE)
png(
  filename = file.path(arg_list$out_dir, "heatmap.png"),
  width = ComplexHeatmap:::width(hm),
  height = ComplexHeatmap:::height(hm),
  units = "mm",
  res = 600
)
draw(hm, annotation_legend_list = lgd_list, merge_legend = TRUE)
dev.off()

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

