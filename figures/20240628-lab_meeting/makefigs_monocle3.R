# ==============================================================================
# These figures were generated for the lab meeting on 2024-06-28. This script
# should be run in the `monocle3` conda environment.
# ==============================================================================

library(tidyverse)
library(monocle3)

# >>> load Monocle object >>>

# import monocle pseudotime results
human_cds <- qs::qread("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240611/human_pseudotime/post_pseudotime_cds.qs")
mouse_cds <- qs::qread("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240611/mouse_pseudotime/post_pseudotime_cds.qs")

# <<< load Monocle object <<<

# set default `cell_size` and `cell_stroke` when plotting to get rid of the ugly
# black border that Monocle makes
cell_size <- 0
cell_stroke <- 0.6

# ------------------------------------------------------------------------------
# human trajectory and pseudotime

# plot cell types
.plt <- plot_cells(
  human_cds,
  reduction = "UMAP",
  color_cells_by = "cell_type_annot",
  group_cells_by = "cell_type_annot",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_principal_points = FALSE
)
ggsave(
  filename = "human_trajectory_cell_type.png",
  plot = .plt,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)

# plot pseudotime
.plt <- plot_cells(
  human_cds,
  reduction = "UMAP",
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_principal_points = TRUE
)
ggsave(
  filename = "human_pseudotime.png",
  plot = .plt,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)

# plot actual age
.plt <- plot_cells(
  human_cds,
  reduction = "UMAP",
  color_cells_by = "human_age",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  label_principal_points = FALSE
)
ggsave(
  filename = "human_age.png",
  plot = .plt,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# mouse trajectory and pseudotime

# plot cell types
.plt <- plot_cells(
  mouse_cds,
  reduction = "UMAP",
  color_cells_by = "cell_type_annot",
  group_cells_by = "cell_type_annot",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_principal_points = FALSE
)
ggsave(
  filename = "mouse_trajectory_cell_type.png",
  plot = .plt,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)

# plot pseudotime
.plt <- plot_cells(
  mouse_cds,
  reduction = "UMAP",
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_principal_points = TRUE
)
ggsave(
  filename = "mouse_pseudotime.png",
  plot = .plt,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)

# plot actual age
.plt <- plot_cells(
  mouse_cds,
  reduction = "UMAP",
  color_cells_by = "mouse_age",
  label_cell_groups = FALSE,
  cell_size = cell_size,
  cell_stroke = cell_stroke,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  label_principal_points = FALSE
)
ggsave(
  filename = "mouse_age.png",
  plot = .plt,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)
