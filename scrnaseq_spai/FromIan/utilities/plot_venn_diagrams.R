#' Plot Venn diagram from two gene lists (or any two sets).
#'
#' @param venn_data A VennPlotData object generated using `RVenn::Venn()` and
#'   `ggVennDiagram::process_data()`.
#' @param label_intersect If true (default), then common elements will be
#'   labelled.
#'
#' @return A ggplot object.
#' 
venn_plot <- function(
  venn_data,
  label_intersect = TRUE
) {
  # make the plot
  .plt <- ggplot() + 
    
    # the set areas
    geom_sf(aes(fill = count), data = venn_region(venn_data)) + 
    
    # the set edges
    geom_sf(color = "black", size = 2, data = venn_setedge(venn_data), show.legend = FALSE) + 
    
    # name of the sets
    geom_sf_text(
      aes(label = name),
      data = venn_setlabel(venn_data),
      size = 6,
      fontface = "bold"
    ) +
    
    # counts
    geom_sf_text(aes(label = count), data = venn_region(venn_data), size = 5)
  
  if (label_intersect) {
    # get elements in the intersection of two sets
    common_elems <- venn_region(venn_data)$item[venn_region(venn_data)$id == "12"] %>% 
      unlist()
    
    # label intersecting genes
    .plt <- .plt + 
      geom_segment(
        aes(x = 500, y = 400, xend = 500, yend = 200),
        color = "red",
        linewidth = 1
      ) + 
      geom_text(
        aes(x = 500, y = 180, label = paste(common_elems, collapse = ", ")),
        data = venn_region(venn_data),
        size = 5,
        vjust = 1
      )
  }
  
  return(.plt)
}


#' Create a `VennPlotData` object from a list of vectors
#'
#' @param sets A list of at least two vectors, each containing a set of elements
#'   (e.g., genes).
#'
#' @return A `VennPlotData` object that can then be passed into `ggVennDiagram`
#'   and `ggplot2` for plotting.
#' 
prep_venn_data <- function(
  sets
) {
  venn_data <- sets %>% 
    Venn() %>% 
    process_data()
  
  return(venn_data)
}
