# ==============================================================================
# Map clusters from Aldinger/Sepp human UBCs to human+mouse UBC integration.
# ==============================================================================

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(ComplexHeatmap)
library(ggalluvial)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # file to Aldinger/Sepp human UBCs and name of metadata column with cluster
  # info to project onto query dataset
  "--ubc_ref",
  nargs = 2,
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # file to integrated human/mouse UBCs and name of metadata column with
  # clusters
  "--ubc_query",
  nargs = 2,
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  # for testing and troubleshooting
  arg_list <- parser$parse_args(c(
    "--ubc_ref", "/.mounts/labs/pailab/public/HumanMouseUBC/data/UBC.Harmony.RDS", "SCT_snn_res.0.5",
    "--ubc_query", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240825/ubc_subset.qs", "subclusters",
    "--out_dir", "/.mounts/labs/pailab/private/projects/HumanMouseUBC/integrated_human_ubc/20250414/"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", arg_list$out_dir))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# functions, set variables

source("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/software/utilities/plotting.R")

# ------------------------------------------------------------------------------
# load seurat objects

srat <- map(
  .x = c(arg_list$ubc_ref[[1]], arg_list$ubc_query[[1]]),
  .f = function(x) {
    ext <- tools::file_ext(x)
    if (grepl("^rds$", ext, ignore.case = TRUE)) {
      message(sprintf("***Reading in %s***", x))
      srat <- readRDS(x)
    } else if (grepl("^qs$", ext, ignore.case = TRUE)) {
      message(sprintf("***Reading in %s***", x))
      srat <- qs::qread(x)
    } else {
      stop(sprintf(
        "***Cannot read %s; please provide a `.rds` or `.qs` file.***",
        x
      ))
    }
    return(srat)
  }
) %>% 
  set_names(c("ubc_ref", "ubc_query"))

# ------------------------------------------------------------------------------
# transfer labels from human-only UBCs to integrated human/mouse UBCs

# find set of anchors between reference and query datasets
message("***Running `FindTransferAnchors`***")
t0 <- Sys.time()
anchors <- FindTransferAnchors(
  reference = srat$ubc_ref,
  query = srat$ubc_query,
  normalization.method = "SCT",
  reference.assay = "SCT",
  query.assay = "SCT",
  dims = 1:25 # default is 1:30
)
t1 <- Sys.time()
message(sprintf(
  "***`FindTransferAnchors` completed in %s***",
  format(t1 - t0, digits = 3)
))

# project PCA from reference to query dataset and predict labels (i.e. clusters)
# in query dataset
message("***Running `TransferData`***")
t0 <- Sys.time()
predictions <- TransferData(
  anchorset = anchors,
  refdata = arg_list$ubc_ref[[2]],
  reference = srat$ubc_ref
)
t1 <- Sys.time()
message(sprintf(
  "***`TransferData` completed in %s***",
  format(t1 - t0, digits = 3)
))

# save table of prediction results
message(sprintf(
  "***Saving predicted labels to %s***",
  file.path(arg_list$out_dir, "predicted_labels.tsv")
))
write.table(
  x = predictions,
  file = file.path(arg_list$out_dir, "predicted_labels.txt"),
  quote = FALSE,
  sep = "\t"
)

# add predicted clusters to reference and query datasets
srat$ubc_ref$predicted_clusters <- predictions[rownames(srat$ubc_ref[[]]), "predicted.id"]
srat$ubc_query$predicted_clusters <- predictions[rownames(srat$ubc_query[[]]), "predicted.id"]

# ------------------------------------------------------------------------------
# make UMAPs

# UMAPs
.plt <- map(
  .x = c("ubc_ref", "ubc_query"),
  .f = \(x) {
    title <- x

    # plot original clustering and predicted clusters
    .plt <- map(
      .x = c(arg_list[[x]][[2]], "predicted_clusters"),
      .f = \(group_by) {
        title <- paste(x, group_by)
        .plt <- DimPlot(
          srat[[x]],
          reduction = "umap",
          group.by = group_by,
          label = TRUE,
          repel = TRUE
        ) + 
          NoLegend() + 
          labs(title = title)
        return(.plt)
      }
    )

    return(.plt)
  }
) %>%
  list_flatten()
.plt <- wrap_plots(.plt, ncol = 2)
ggsave(
  plot = .plt,
  filename = "umap_clusters.png",
  path = arg_list$out_dir,
  width = 10,
  height = 10,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# bar plot of predicted clusters

p1 <- cluster_barplot(
  object = srat$ubc_query,
  split.by = "species",
  group.by = "predicted_clusters",
  position = "stack"
)
p2 <- cluster_barplot(
  object = srat$ubc_query,
  split.by = "predicted_clusters",
  group.by = "species",
  position = "stack"
)
layout <- "AAAB"
.plt <- wrap_plots(p1, p2, design = layout)
ggsave(
  plot = .plt,
  filename = "barplot_pred_species.png",
  path = arg_list$out_dir,
  width = 10,
  height = 4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# compare original and predicted clusters in reference dataset (to see how well
# TransferData is performing)

message(sprintf(
  "***%s/%s cells in reference dataset were not in query dataset***",
  sum(is.na(srat$ubc_ref$predicted_clusters)),
  nrow(srat$ubc_ref[[]])
))

# get original clusters and predicted clusters from reference
ref_clust_df <- srat$ubc_ref[[]] %>%
  select(
    cluster = arg_list$ubc_ref[[2]],
    predicted_clusters
  ) %>%
  # filter out cells that were not in query dataset (i.e. no predicted cluster)
  filter(!is.na(predicted_clusters))

# matrix of pairwise Jaccard similarity (see
# https://bioconductor.org/books/3.14/OSCA.advanced/clustering-redux.html)
mat <- bluster::linkClustersMatrix(
  x = ref_clust_df$cluster,
  y = ref_clust_df$predicted_clusters,
  denominator = "union"
)
rownames(mat) <- paste0("orig_", rownames(mat))
colnames(mat) <- paste0("pred_", colnames(mat))

# colour palette
col_fun <- circlize::colorRamp2(
  breaks = seq(0, 1, by = 0.05),
  colors = scales::pal_viridis()(21)
)
# show Jaccard similarity in heatmap
cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 8))
}

# make heatmap and save
hm <- Heatmap(
  mat,
  name = "Jaccard\nsimilarity",
  col = col_fun,
  cell_fun = cell_fun,
  row_title = "original clusters",
  row_title_side = "left",
  row_names_side = "left",
  column_title = "predicted clusters",
  column_title_side = "bottom",
  column_names_side = "bottom",
  column_names_rot = 45,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  width = ncol(mat) * unit(0.5, "in"),
  height = nrow(mat) * unit(0.5, "in"),
)
hm <- draw( # draw to get dimensions
  hm,
  column_title = "Human-only UBCs (Aldinger/Sepp)",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  padding = unit(c(5.5, 5.5, 5.5, 11), "pt")
)
png(
  filename = file.path(arg_list$out_dir, "ref_dataset_jaccard_heatmap.png"),
  width = ComplexHeatmap:::width(hm),
  height = ComplexHeatmap:::height(hm),
  units = "mm",
  res = 600
)
draw(hm)
dev.off()

# save matrix
write.table(
  x = mat,
  file = file.path(arg_list$out_dir, "ref_dataset_jaccard_matrix.txt"),
  quote = FALSE,
  sep = "\t"
)

# ------------------------------------------------------------------------------
# compare original and predicted clusters in query dataset

# get original clusters and predicted clusters from query
query_clust_df <- srat$ubc_query[[]] %>%
  select(
    cluster = arg_list$ubc_query[[2]],
    predicted_clusters,
    species
  )

# matrix of pairwise Jaccard similarity
mat <- bluster::linkClustersMatrix(
  x = query_clust_df$cluster,
  y = query_clust_df$predicted_clusters,
  denominator = "union"
)
rownames(mat) <- paste0("orig_", rownames(mat))
colnames(mat) <- paste0("pred_", colnames(mat))

# show Jaccard similarity in heatmap
cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 8, col = "grey75"))
}

# make heatmap and save
hm <- Heatmap(
  mat,
  name = "Jaccard\nsimilarity",
  col = col_fun,
  cell_fun = cell_fun,
  row_title = "original clusters",
  row_title_side = "left",
  row_names_side = "left",
  column_title = "predicted clusters",
  column_title_side = "bottom",
  column_names_side = "bottom",
  column_names_rot = 45,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  width = ncol(mat) * unit(0.5, "in"),
  height = nrow(mat) * unit(0.5, "in"),
)
hm <- draw( # draw to get dimensions
  hm,
  column_title = "Integrated human/mouse UBCs",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  padding = unit(c(5.5, 5.5, 5.5, 11), "pt")
)
png(
  filename = file.path(arg_list$out_dir, "query_dataset_jaccard_heatmap.png"),
  width = ComplexHeatmap:::width(hm),
  height = ComplexHeatmap:::height(hm),
  units = "mm",
  res = 600
)
draw(hm)
dev.off()

# save matrix
write.table(
  x = mat,
  file = file.path(arg_list$out_dir, "query_dataset_jaccard_matrix.txt"),
  quote = FALSE,
  sep = "\t"
)

# ------------------------------------------------------------------------------
# compare original and predicted clusters in query dataset by species

species_list <- c("human", "mouse")
hm <- list()
show_heatmap_legend <- TRUE

for (species in species_list) {
  mat <- bluster::linkClustersMatrix(
    x = query_clust_df$cluster[query_clust_df$species == species],
    y = query_clust_df$predicted_clusters[query_clust_df$species == species],
    denominator = "union"
  )

  # add back missing columns (if any) to ensure consistent dimensions
  all_clusters <- unique(query_clust_df$predicted_clusters) %>% sort()
  missing_cols <- setdiff(all_clusters, colnames(mat))
  if (length(missing_cols) > 0) {
    mat <- cbind(
      mat,
      matrix(0, nrow = nrow(mat), ncol = length(missing_cols), dimnames = list(NULL, missing_cols))
    )
    mat <- mat[, all_clusters]
  }

  rownames(mat) <- paste0("orig_", rownames(mat))
  colnames(mat) <- paste0("pred_", colnames(mat))

  hm[[species]] <- Heatmap(
    mat,
    name = if (show_heatmap_legend) {"Jaccard\nsimilarity"} else {species},
    col = col_fun,
    cell_fun = cell_fun,
    row_title = sprintf("original clusters\n(%s only)", species),
    row_title_side = "left",
    row_names_side = "left",
    column_names_side = "bottom",
    column_names_rot = 45,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    width = ncol(mat) * unit(0.5, "in"),
    height = nrow(mat) * unit(0.5, "in"),
    show_heatmap_legend = show_heatmap_legend
  )

  show_heatmap_legend <- FALSE

  # save matrix
  write.table(
    x = mat,
    file = file.path(arg_list$out_dir, sprintf("query_dataset_jaccard_matrix_%s.txt", species)),
    quote = FALSE,
    sep = "\t"
  )
}

hm_list <- hm$human %v% hm$mouse

hm_list <- draw( # draw to get dimensions
  hm_list,
  column_title = "Integrated human/mouse UBCs\n(split by species)",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  padding = unit(c(5.5, 5.5, 5.5, 11), "pt")
)
png(
  filename = file.path(arg_list$out_dir, "query_dataset_jaccard_heatmap_species.png"),
  width = ComplexHeatmap:::width(hm_list),
  height = ComplexHeatmap:::height(hm_list),
  units = "mm",
  res = 600
)
draw(hm_list)
dev.off()


# ---
# make alluvial plot with ggalluvial

plt_data <- query_clust_df %>%
  group_by(cluster, predicted_clusters, species) %>%
  summarise(freq = n())

.plt <- ggplot(plt_data, aes(y = freq, axis1 = cluster, axis2 = predicted_clusters)) + 
  geom_alluvium(aes(fill = species), width = 0.25) + 
  geom_stratum(width = 0.25) + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) + 
  scale_x_discrete(
    limits = c("original clusters", "predicted clusters"),
    expand = expansion(mult = 0.25)
  ) + 
  labs(
    title = "Original vs predicted clusters\nin integrated human/mouse UBCs",
    y = "number of cells"
  ) +
  theme_bw()
ggsave(
  plot = .plt,
  filename = "alluvial_original_vs_predicted_clusters.png",
  path = arg_list$out_dir,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())
