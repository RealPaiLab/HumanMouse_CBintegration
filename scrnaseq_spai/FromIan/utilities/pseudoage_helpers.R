# ==============================================================================
# Functions to help calculate pseudoages of each cell to determine the
# correspondence of developmental stages across mouse and human (see methods
# from Sepp et al. 2024).
# ==============================================================================


#' Find k-nearest neighbours of each cell.
#'
#' @description Note that this function uses the Python `pynndescent` package
#'   which only look for **approximate** nearest neighbours, not exact nearest
#'   neighbours.
#'
#' @param query_mtx Cell-by-embedding matrix of the cells for which to find the
#'   nearest neighbours in the reference matrix (`ref_mtx`).
#' @param ref_mtx Cell-by-embedding matrix of the reference cells.
#' @param k Number of nearest neighbour cells to find.
#' @param metric Distance metric to use for determining nearest neighbours
#'   (passed to the `NNDescent` function).
#'
#' @return A list of two matrices: (1) a matrix with the query cells in the rows
#'   and the index of the `k` nearest neighbours in the columns; (2) a matrix
#'   with the query cells in the rows and the distances to the `k` nearest
#'   neighbours in the columns. See the `query` method in the Python
#'   `pynndescent` package for more details.
#'
find_knn <- function(
  query_mtx,
  ref_mtx,
  k = 25,
  metric = "cosine"
) {
  # import the python package pynndescent
  reticulate::use_condaenv("scrnaseq_parallel_env")
  pynndescent <- reticulate::import("pynndescent")

  # build index
  message("***Building nearest neighbour index***")
  index <- pynndescent$NNDescent(ref_mtx, metric = metric)

  # find (approximate) nearest neighbours
  message("***Finding (approximate) nearest neighbours***")
  nn <- index$query(query_mtx, k = as.integer(k))
  names(nn) <- c("nn_index", "nn_dist")

  # the `query` returns zero-based indices (cuz Python);
  # need to convert back to one-based (for R)
  nn$nn_index <- nn$nn_index + 1

  return(nn)
}


#' Calculate pseudostage of each cell in a Seurat object.
#'
#' @param srat Query Seurat object for which to calculate pseudostage.
#' @param ref_srat Reference Seurat object with the age index for each cell.
#' @param ref_stages Name of metadata column in reference Seurat object
#'   containing the age index for each cell.
#' @param reduction Name of the reduction in Seurat object to use for
#'   calculating nearest neighbour distance.
#' @param k Number of nearest neighbour cells to find.
#' @param metric Distance metric to use for determining nearest neighbours
#'   (passed to the `NNDescent` function).
#'
#' @return A named vector containing the pseudostage of each cell.
#'
get_pseudostage <- function(
  srat,
  ref_srat,
  ref_stages,
  reduction,
  k = 25,
  metric = "cosine"
) {
  message(sprintf(
    "***Using %s embeddings in the %s assay***",
    reduction,
    srat@active.assay
  ))

  # get indices (row numbers) of the k nearest neighbours
  nn_index <- find_knn(
    query_mtx = Embeddings(srat, reduction = reduction),
    ref_mtx = Embeddings(ref_srat, reduction = reduction),
    k = k,
    metric = metric
  ) %>%
    pluck("nn_index")

  # get mean age index (stage) of the k nearest neighbours
  message("***Calculating mean pseudoage of each cell***")
  pseudostage <- apply(
    X = nn_index,
    MARGIN = 1, # 1 for row of matrix (i.e. per cell)
    FUN = function(X) {
      # get column containing stages and convert to vector
      out <- ref_srat[[ref_stages]][[1]] %>%
        # subset stages of nearest neighbours
        { .[X] } %>%
        # calculate mean
        mean()

      return(out)
    }
  )
  names(pseudostage) <- rownames(srat[[]])

  # return output
  return(pseudostage)
}


#' Bin cells by pseudoage and calculate proportion of cells in each bin at each
#' developmental timepoint. See notes from 2024-05-07 and 2024-05-15 for more
#' details about this function.
#'
#' @param md Metadata table from Seurat containing the `pseudoage` and the
#'   actual human/mouse ages of each cell.
#' @param bin_col Name of the pseudoage column in the metadata table to use
#'   for binning the pseudoages. Defaults to "pseudoage".
#'
#' @return A dataframe with the proportion of cells in each pseudoage bin
#'   (columns) for each developmental timepoint (rows). Rows should sum to 1.
#'
bin_pseudoage <- function(
  md,
  bin_col = "pseudoage"
) {
  # bin cells by pseudoage
  x <- mutate(
    md,
    pseudoage_bin = cut(.data[[bin_col]], breaks = 50, labels = FALSE) %>% as.factor()
  ) %>%
    # keep only certain columns
    select(
      species,
      mouse_age,
      human_age,
      pseudoage_bin
    )

  # calculate proportion of cells in each bin for each developmental timepoint
  x <- x %>%
    mutate(
      all_age = case_when(
        species == "mouse" ~ mouse_age,
        species == "human" ~ human_age
      )
    ) %>%
    group_by(species, all_age, pseudoage_bin) %>%
    count() %>%
    ungroup(pseudoage_bin) %>%
    mutate(
      prop = n / sum(n)
    )

  # reformat so each developmental timepoint is in a row and each bin is in a
  # column
  x <- x %>%
    select(-n) %>%
    pivot_wider(
      names_from = pseudoage_bin,
      values_from = prop,
      names_sort = TRUE,
      values_fill = 0
    ) %>%
    ungroup() %>%
    select(-species) %>%
    column_to_rownames("all_age") %>%

  return(x)
}
