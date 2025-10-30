#' Generate Seurat FeaturePlot and save as PDF.
#' 
#' @param object Seurat object
#' @param features Vector of features to plot
#' @param plot_args List of arguments to be passed to `FeaturePlot`
#' @param save_args List of arguments to be passed to `ggsave`
#' 
save_feature_plots <- function(
  object,
  features,
  plot_args = list(),
  save_args = list()
) {
  # loop through features
  for (feature in features) {
    
    # make plot
    feature_plot <- do.call(
      FeaturePlot,
      args = c(list(object = object, features = feature), plot_args)
    )
    
    # save plot
    filename <- paste0(feature, ".pdf")
    do.call(
      ggsave,
      args = c(list(filename = filename, plot = feature_plot), save_args)
    )
  }
}

#' Generate VlnPlot and save as PDF.
#' 
#' @param object Seurat object
#' @param features Vector of features to plot
#' @param plot_args List of arguments to be passed to `VlnPlot`
#' @param save_args List of arguments to be passed to `ggsave`
#' 
save_vln_plots <- function(
  object,
  features,
  plot_args = list(),
  save_args = list()
) {
  # loop through features
  for (feature in features) {
    
    # make plot
    feature_plot <- do.call(
      VlnPlot,
      args = c(list(object = object, features = feature), plot_args)
    )
    
    # save plot
    filename <- paste0(feature, ".pdf")
    do.call(
      ggsave,
      args = c(list(filename = filename, plot = feature_plot), save_args)
    )
  }
}
