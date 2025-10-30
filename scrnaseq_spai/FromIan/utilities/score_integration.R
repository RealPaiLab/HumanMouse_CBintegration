# ==============================================================================
# Functions to score the integration of cells from different species by
# comparing labelled cell types.
# ==============================================================================


#' Calculate the distance between species centroid embeddings split by cell
#' type.
#' 
#' @param embeddings
#' @param metadata
#' @param delta Metadata on which to compute the absolute difference in centroid
#'   embeddings. Default is "species" which returns the absolute difference
#'   between the centroids of the human cells and the mouse cells.
#' @param by Metadata to group the delta centroids by. Default is the
#'   "common_cell_type" which calculates the absolute delta centroid for each
#'   cell type.
#' @param filter_out Cell types (or other group) to filter out.
#' @param integ_method String specifying the integration method that was used.
#'   If NA, no column is added to the output data frame; otherwise, a column
#'   with the string is added.
#'   
#' @return Data frame with the distance between human and mouse centroids (delta
#'   centroids) grouped by the `by` parameter. By default, returns the delta
#'   centroids for each common cell type.
#'
delta_centroids <- function(
  embeddings,
  metadata,
  # delta = "species",
  by = "common_cell_type",
  filter_out = c("other/missing"), 
  integ_method = NA
) {
  # get human centroids
  human <- calc_centroids(
    embeddings = embeddings[metadata["species"] == "human", ],
    metadata = metadata[metadata["species"] == "human", ],
    by = by,
    filter_out = filter_out
  )
  
  # get mouse centroids
  mouse <- calc_centroids(
    embeddings = embeddings[metadata["species"] == "mouse", ],
    metadata = metadata[metadata["species"] == "mouse", ],
    by = by,
    filter_out = filter_out
  )
  
  # make sure cell types (or other grouping) are the same between human and
  # mouse (e.g. if RL is present in human but not in mouse, RL needs to be
  # filtered out)
  if (!is_empty(setdiff(colnames(human), colnames(mouse)))) {
    stop(sprintf("`%s` in human and mouse do not match. Try using `filter_out` to remove groups that are missing from one species.", by))
  }
  
  # calculate distance between human and mouse centroid for each of the common cell types
  all_dist <- data.frame()
  for (cname in colnames(human)) {
    # subset by cell type (or other grouping) and calculate Euclidian distance
    delta <- rbind(human[, cname], mouse[, cname]) %>% 
      dist(.)
    
    # save to data frame
    all_dist <- rbind(all_dist, c(cname, delta))
  }
  
  # set column names of the data frame containing the distances between centroids
  colnames(all_dist) <- c(by, "delta")
  
  # add integration method as column
  if (!is.na(integ_method)) {
    all_dist <- mutate(
      all_dist,
      integ_method = integ_method,
      .before = 1
    )
  }
  
  # convert distance column to numeric (not sure why it ends up as a character)
  all_dist <- mutate(
    all_dist,
    delta = as.numeric(delta)
  )
  
  return(all_dist)
}



#' Calculate the centroid coordinates grouped by `by`.
#'
#' @param embeddings
#' @param metadata
#' @param by
#' @param filter_out Cell types (or other group) to filter out.
#'
#' @return A data frame of centroids with n rows and g columns, where g is the
#'   number of groups and n is the n-dimensional coordinates of the centroid.
#'   
calc_centroids <- function(
  embeddings,
  metadata,
  by = "common_cell_type",
  filter_out = c("other/missing")
) {
  centroids <- unique(metadata[[by]]) %>% 
    # calculate the centroid (`colMeans`) for each cell type in metadata[[by]]
    purrr::map_dfr(~ colMeans(embeddings[metadata[by] == ., ])) %>% 
    t(.) %>% 
    as.data.frame(.) %>% 
    # set the cell types (or other grouping) as the column names
    `colnames<-`(unique(metadata[[by]])) %>% 
    # filter out specific cell types
    select(!any_of(filter_out))
  
  return(centroids)
}
