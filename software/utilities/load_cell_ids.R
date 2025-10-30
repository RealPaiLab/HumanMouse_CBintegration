#' Load cell IDs for previously annotated cell types.
#' 
#' Use `TRUE`/`FALSE` to indicate whether or not to return the cell IDs for each
#' cell type.
#'
#' @param rl_vz Cell IDs for the Hendrikse-labelled RL-VZ cells (human cells only).
#' @param rl_svz Cell IDs for the Hendrikse-labelled RL-SVZ cells (human cells only).
#' @param integ_c7 Cell IDs for cluster 7 ("common UBCs") in the inital
#'   human/mouse integration (human and mouse cells).
#' @param integ_c19 Cell IDs for cluster 19 ("human-specific UBCs") in the inital
#'   human/mouse integration (human and a few mouse cells).
#' @param integ_c20 Cell IDs for cluster 20 ("human-specific UBCs") in the inital
#'   human/mouse integration (human and a few mouse cells).
#'
#' @return List with the cell IDs for all requested cell types.
#'
load_cell_ids <- function(
  rl_vz = TRUE,
  rl_svz = TRUE,
  integ_c7 = TRUE,
  integ_c19 = TRUE,
  integ_c20 = TRUE
) {
  main_dir <- "/.mounts/labs/pailab/private/llau/results/integrated"
  cell_ids <- list()

  get_cell_ids <- function(f) {
    ids <- read.csv(file = f, header = TRUE, row.names = 1) %>%
      dplyr::pull(x)
    return(ids)
  }

  if (rl_vz) {
    cell_ids$rl_vz <- file.path(main_dir, "20240509/cca/vz_cells.csv") %>%
      get_cell_ids()
  }

  if (rl_svz) {
    cell_ids$rl_svz <- file.path(main_dir, "20240509/cca/svz_cells.csv") %>%
      get_cell_ids()
  }

  if (integ_c7) {
    cell_ids$cluster_7 <- file.path(main_dir, "20240516/fc_subset_analysis/cluster_7_cells.csv") %>%
      get_cell_ids()
  }

  if (integ_c19) {
    cell_ids$cluster_19 <- file.path(main_dir, "20240516/fc_subset_analysis/cluster_19_cells.csv") %>%
      get_cell_ids()
  }

  if (integ_c20) {
    cell_ids$cluster_20 <- file.path(main_dir, "20240516/fc_subset_analysis/cluster_20_cells.csv") %>%
      get_cell_ids()
  }

  return(cell_ids)
}
