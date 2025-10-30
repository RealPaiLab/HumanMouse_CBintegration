#' @import Seurat
#' @import qs
#' @import tidyverse

#' Creates jaccard similarity matrix between clusters of two integrations
#' 
#' @param srat_1 First Seurat object.
#' @param res_1 First Seurat object's resolution.
#' @param srat_2 Second Seurat object.
#' @param res_2 Second Seurat object's resolution.

jaccard_matrix <- function(
    srat_1,
    res_1,
    srat_2,
    res_2
){
  res_1_col <- paste0("snn_res.", res_1)
  res_2_col <- paste0("snn_res.", res_2)

  x <- list()
  for(level in as.integer(levels(srat_1@meta.data[[res_1_col]]))) {
    cells <- rownames(srat_1@meta.data[srat_1@meta.data[[res_1_col]] == level, ])
    x[[length(x) + 1]] <- cells
  }

  y <- list()
  for(level in as.integer(levels(srat_2@meta.data[[res_2_col]]))) {
    cells <- rownames(srat_2@meta.data[srat_2@meta.data[[res_2_col]] == level, ])
    y[[length(y) + 1]] <- cells
  }

  jaccard_matrix <- outer(x, y, FUN = Vectorize(jaccard_sim))
  rownames(jaccard_matrix) <- paste0("UBC_", as.integer(levels(srat_1@meta.data[[res_1_col]])))
  colnames(jaccard_matrix) <- paste0("UBC_", as.integer(levels(srat_2@meta.data[[res_2_col]])))

  return(jaccard_matrix)
}

jaccard_sim <- function(
    cells_a,
    cells_b
){
  intersection <- length(intersect(cells_a, cells_b))
  union <- length(union(cells_a, cells_b))
  jaccard <- intersection / union
  return(jaccard)
}