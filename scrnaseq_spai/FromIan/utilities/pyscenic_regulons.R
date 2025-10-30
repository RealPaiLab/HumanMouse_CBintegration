#' Get top scoring regulons by RSS for the specified cell type.
#'
#' @param cell_type Cell type to use when sorting the RSS
#' @param rss Dataframe of regulon specificity scores with a cell type in each
#'   row and a TF/regulator in each column.
#' @param num_regulons Number of regulons to return. Defaults to the top 50%.
#'   `Inf` returns all regulons.
#'
#' @return Vector of regulon names (i.e., regulators) from highest to lowest RSS
#'   for the cell type.
#' 
get_top_regulons <- function(
  cell_type,
  rss,
  num_regulons = ceiling(ncol(rss)/2)
) {
  scores <- data.matrix(rss)[cell_type, ]
  
  # sort scores then get the TF associated with the score
  top_regs <- scores %>% 
    sort(decreasing = TRUE) %>% 
    names() %>% 
    head(n = num_regulons)
  
  return(top_regs)
}
