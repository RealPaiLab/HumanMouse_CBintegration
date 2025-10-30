#' @import Seurat
#' @import tidyverse

#' Integration of list of datasets using RPCA method
#' 
#' @param integ_list List of Seurat dataset objects in first column ($dataset_list), and integration featurtes in second column (integ_feat)
#' @param resolutions Vector of numerical characters that specify the clustering resolutions. Default: (0.1, 0.3, 0.5, 0.7, 0.9)
#' @param do_tsne Boolean indicating to do T-SNE plots
#' @param out_directory String of output directory
#' @param log_file String of logfile location 
#' @param use_rds Boolean indicating to write .rds files


rpca_integration <- function(
  integ_list,
  resolutions =  seq(0.1, 0.9, by = 0.2),
  do_tsne = FALSE,
  out_directory,
  log_file = NULL,
  use_rds = FALSE
) {
  # If want to integrate in series, all outputs and errors should print to main log
  if(!is.null(log_file)){
    sink(log_file, type = "output")
    sink(log_file, type = "message", append = TRUE)
  }

  plan(multicore)

  # Extracting necessary integration variables
  dataset_list <- integ_list$dataset_list
  features <- integ_list$integ_feat

  # Keeping track of k.anchor value
  k_anchor = 60
  cat("K.anchor value:", k_anchor, "\n")

  # find anchors and integrate
  print("Finding RPCA integration anchors")
  anchors_rpca <- FindIntegrationAnchors(dataset_list, normalization.method = "SCT", anchor.features = features, reduction = "rpca", k.anchor = k_anchor) 
  print("Performing RPCA integration")
  integ_srat_rpca <- IntegrateData(anchorset = anchors_rpca, normalization.method = "SCT")

  # run PCA, check PCs
  print("Running PCA for RPCA integration")
  integ_srat_rpca <- RunPCA(integ_srat_rpca, features = VariableFeatures(integ_srat_rpca), npcs = 100)
  rpca_dim_heatmap <- DimHeatmap(integ_srat_rpca, dims = 1:15, cells = 500, fast = FALSE)
  ggsave(filename = "RPCA_dim_heatmap.png", plot = rpca_dim_heatmap, path = out_directory, 
          width = 10, height = 15, units = "in", dpi = 600)
  rpca_elbow_plot <- ElbowPlot(integ_srat_rpca, ndims = 100) + theme_classic()
  ggsave(filename = "RPCA_elbow_plot.png", plot = rpca_elbow_plot, path = out_directory, 
          width = 10, height = 7.5, units = "in", dpi = 600)

  # set ndims for clustering and dimensional reductions (based on elbow plot)
  ndims <- 50

  # Run TSNE and UMAP
  if(do_tsne){
    print("Running TSNE for RPCA integration")
    integ_srat_rpca <- RunTSNE(integ_srat_rpca, reduction = "pca", dims = 1:ndims, num_threads = 32)
  }
  print("Running UMAP for RPCA integration")
  integ_srat_rpca <- RunUMAP(integ_srat_rpca, reduction = "pca", dims = 1:ndims)

  # cluster cells
  print("Clustering RPCA integration")
  source("scrnaseq_Leo/utilities/cluster_multiple_res.R")
  integ_srat_rpca <- FindNeighbors(integ_srat_rpca, dims = 1:ndims) %>%
    cluster_multiple_res(resolutions, "RPCA", out_directory)

  # Normalizing data slot in RNA assay. 
  # Some datasets don't have their data slot in RNA normalized, causing dowstream visualization issues
  integ_srat_rpca <- NormalizeData(integ_srat_rpca, assay = "RNA")

  if(use_rds){
    # Save integrated dataset as a .rds file
    print("Saving RPCA integration as .rds file")
    file_name = paste(format(Sys.Date(), "%Y%m%d"), "rpca_integ.rds", sep = "_")
    saveRDS(object = integ_srat_rpca, file = file.path(out_directory, file_name))
  } else {
    # Save integrated dataset as a .qs file
    print("Saving RPCA integration as .qs file")
    file_name = paste(format(Sys.Date(), "%Y%m%d"), "rpca_integ.qs", sep = "_")
    qsave(integ_srat_rpca, file.path(out_directory, file_name))
  }
  
  message("Dataset successfully underwent RPCA integration. The integrated dataset is called: ", file_name)

  source("scrnaseq_Leo/pipeline_scripts/visualization.R")
  visualization(integ_srat_rpca, ndims = ndims, resolutions = resolutions, integ_type = "RPCA", do_tsne = do_tsne, out_directory = out_directory)

  if(!is.null(log_file)){
    sink(type = "output")
    sink(type = "message")
  }
}