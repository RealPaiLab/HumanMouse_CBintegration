#' @import Seurat
#' @import tidyverse

#' Integration of list of datasets using CCA method
#' 
#' @param integ_list List of Seurat dataset objects in first column ($dataset_list), and integration featurtes in second column (integ_feat)
#' @param resolutions Vector of numerical characters that specify the clustering resolutions. Default: (0.1, 0.3, 0.5, 0.7, 0.9)
#' @param do_tsne Boolean indicating to do T-SNE plots
#' @param out_directory String of output directory
#' @param log_file String of logfile location 
#' @param use_rds Boolean indicating to write .rds files


cca_integration <- function(
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

  # find anchors and integrate
  print("Finding CCA integration anchors")
  anchors_cca <- FindIntegrationAnchors(dataset_list, normalization.method = "SCT", anchor.features = features, reduction = "cca")
  print("Performing CCA integration")
  integ_srat_cca <- IntegrateData(anchorset = anchors_cca, normalization.method = "SCT")

  # run PCA, check PCs
  print("Running PCA for CCA integration")
  integ_srat_cca <- RunPCA(integ_srat_cca, features = VariableFeatures(integ_srat_cca), npcs = 100)
  cca_dim_heatmap <- DimHeatmap(integ_srat_cca, dims = 1:15, cells = 500, fast = FALSE)
  ggsave(filename = "CCA_dim_heatmap.png", plot = cca_dim_heatmap, path = out_directory, 
          width = 10, height = 15, units = "in", dpi = 600)
  cca_elbow_plot <- ElbowPlot(integ_srat_cca, ndims = 100) + theme_classic()
  ggsave(filename = "CCA_elbow_plot.png", plot = cca_elbow_plot, path = out_directory, 
          width = 10, height = 7.5, units = "in", dpi = 600)

  # set ndims for clustering and dimensional reductions (based on elbow plot)
  ndims <- 50

  # Run TSNE and UMAP
  if(do_tsne){
    print("Running TSNE for CCA integration")
    integ_srat_cca <- RunTSNE(integ_srat_cca, reduction = "pca", dims = 1:ndims, num_threads = 32)
  }
  print("Running UMAP for CCA integration")
  integ_srat_cca <- RunUMAP(integ_srat_cca, reduction = "pca", dims = 1:ndims)

  # cluster cells
  print("Clustering CCA integration")
  source("scrnaseq_Leo/utilities/cluster_multiple_res.R")
  integ_srat_cca <- FindNeighbors(integ_srat_cca, dims = 1:ndims) %>%
    cluster_multiple_res(resolutions, 'CCA', out_directory)
  
  # Normalizing data slot in RNA assay. 
  # Some datasets don't have their data slot in RNA normalized, causing dowstream visualization issues
  integ_srat_cca <- NormalizeData(integ_srat_cca, assay = "RNA")

  if(use_rds){
    # Save integrated dataset as a .rds file
    print("Saving CCA integration as .rds file")
    file_name = paste(format(Sys.Date(), "%Y%m%d"), "cca_integ.rds", sep = "_")
    saveRDS(object = integ_srat_cca, file = file.path(out_directory, file_name))
  } else {
    # Save integrated dataset as a .qs file
    print("Saving CCA integration as .qs file")
    file_name = paste(format(Sys.Date(), "%Y%m%d"), "cca_integ.qs", sep = "_")
    qsave(integ_srat_cca, file.path(out_directory, file_name))
  }

  message("Dataset successfully underwent CCA integration. The integrated dataset is called: ", file_name)
  
  source("scrnaseq_Leo/pipeline_scripts/visualization.R")
  visualization(integ_srat_cca, ndims = ndims, resolutions = resolutions, integ_type = "CCA", do_tsne = do_tsne, out_directory = out_directory)

  if(!is.null(log_file)){
    sink(type = "output")
    sink(type = "message")
  }
}