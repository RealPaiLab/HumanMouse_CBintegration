#' @import ggplot2

#' Stacked bar plot showing composition of each cluster.
#' 
#' @param object Seurat object
#' @param split.by metadata variable to split each bar by
#' @param group.by metadata variable to group each bar by
#' @param position parameter passed to `geom_col`, either "stack" or "fill"
#' @param width width of bars, passed to `geom_col`
#' @param filter_data string to filter out some `split.by` and/or `group.by`
#'   variables from the plot, e.g. `filter_data = "seurat_clusters %in% c(1, 5)"`
#'   will only show the bar plot for seurat_clusters 1 and 5
#' 
cluster_barplot <- function(
  object,
  split.by,
  group.by = "seurat_clusters",
  position = "stack",
  width = 0.75,
  filter_data = NULL
) {
  # deparse(substitute(foo)) takes foo and converts it to a string
  num_cells <- table(
    object@meta.data[[split.by]],
    object@meta.data[[group.by]]
  ) %>%
    as.data.frame(.)
  
  # rename columns
  colnames(num_cells) <- c(split.by, group.by, "freq")
  
  # filter in/out specific rows
  if (!is.null(filter_data)) {
    num_cells <- dplyr::filter(
      num_cells,
      !! rlang::parse_expr(filter_data)
    )
  }
  
  plt <- ggplot(
    data = num_cells,
    aes(x = .data[[group.by]], y = freq, fill = .data[[split.by]])
  ) +
    geom_col(position = position, width = width) + 
    labs(y = if (position == "stack") {
      "number of cells"
    } else if (position == "fill") {
      "proportion"
    } else {
      position
    }) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    theme_light() + 
    theme(axis.text = element_text(colour = "black"))
  
  return(plt)
}
