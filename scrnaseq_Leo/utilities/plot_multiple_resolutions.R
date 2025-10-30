#' @import Seurat
#' @import tidyverse

#' Run multiple clustering resolutions on raw data after SCTransform() and before RunPCA() using FindCluster().
#' 
#' @param object Seurat object
#' @param reduction dimensionality reduction method
#' @param dims number of dimensions to use
#' @param resolutions resolutions at which to cluster at (in a vector)
#' @param label include label on plot


plot_multiple_resolutions <- function(
  object,
  reduction = "umap",
  dims = 1:30,
  resolutions,
  label = TRUE 
) {
    #Processing
    data_srat <- RunPCA(object) %>%
        FindNeighbors(dims = dims)
    
    plot_list <- list()
    
    for (i in resolutions) {
        message("Generating clustering at ", i, " resolution.")
        data_srat <- FindClusters(data_srat, resolution = i) 
    }
    
    data_srat <- RunUMAP(data_srat, dims = dims)
    c <- 1

    for (i in resolutions) {
        pattern <- paste0(".*_snn_res.", i)
        matching_column <- grep(pattern, colnames(data_srat@meta.data), value = TRUE)
        plt <- DimPlot(data_srat, reduction = "umap", group.by = matching_column, label = label, label.size = 3) + ggtitle(paste("Resolution:", i))
        plot_list[[c]] <- plt
        c <- c + 1
    }
  
  return(plot_list)
}
