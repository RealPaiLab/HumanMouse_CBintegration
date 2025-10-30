#' @import Seurat
#' @import clustree

#' Run multiple clustering resolutions on raw data after SCTransform() and RunPCA() using FindCluster().
#' 
#' @param dataset_srat Seurat object containing dataset
#' @param resolutions Resolutions at which to cluster at (in a vector)
#' @param integ_type String indicating integration type (CCA, Harmony, RPCA)
#' @param out_directory String that specifies output directory


cluster_multiple_res <- function(
  dataset_srat,
  resolutions,
  integ_type,
  out_directory
) {
    # Perform clustering at multiple resolutions
    for (i in resolutions) {
        message("Generating clustering at ", i, " resolution.")
        dataset_srat <- FindClusters(dataset_srat, resolution = i)
        pattern <- paste0(".*_snn_res.", i)
        matching_column <- grep(pattern, colnames(dataset_srat@meta.data), value = TRUE)
        new_name <- paste0("snn_res.", i)
        dataset_srat@meta.data[[new_name]] <- dataset_srat@meta.data[[matching_column]]
        dataset_srat@meta.data[[matching_column]] <- NULL
    }

    # Extract and visualize the clusters
    #tree <- clustree(dataset_srat, prefix = "snn_res.")
    #tree_plot_name = paste0(integ_type, "_tree_plot.png")
    #ggsave(filename = tree_plot_name, plot = tree, path = out_directory, 
            #width = 15, height = 15, units = "in", dpi = 600)

    return(dataset_srat)
}
