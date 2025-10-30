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


harmony_integration <- function(
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

  # Extracting necessary integration variables
  dataset_list <- integ_list$dataset_list
  features <- integ_list$integ_feat
  
  print("Merging datasets")
  merge_srat <- Reduce(merge, dataset_list)

  print("Setting variable features to common genes")
  # set variable features to common genes
  VariableFeatures(merge_srat) <- features

  # run PCA
  print("Running PCA for Harmony integration")
  merge_srat <- RunPCA(merge_srat, features = VariableFeatures(merge_srat), npcs = 100)

  # run harmony
  print("Running Harmony integration")

  # Keeping track of theta value
  theta = c(2, 5)
  cat("Theta value:", theta, "\n")

  integ_srat_harm <- RunHarmony(
    object = merge_srat,
    group.by.vars = c("species", "dataset_name"),
    assay.use = "SCT",
    project.dim = FALSE,
    theta = theta
  )

  harmony_elbow_plot <- ElbowPlot(integ_srat_harm, ndims = 50, reduction = "harmony") + theme_classic()
  ggsave(filename = "Harmony_elbow_plot.png", plot = harmony_elbow_plot, path = out_directory, 
          width = 10, height = 7.5, units = "in", dpi = 600)

  # set ndims for clustering and dimensional reductions (based on elbow plot)
  ndims <- 50

  # Run TSNE and UMAP
  if(do_tsne){
    print("Running TSNE for Harmony integration")
    integ_srat_harm <- RunTSNE(integ_srat_harm, reduction = "harmony", dims = 1:ndims, num_threads = 32)
  }
  print("Running UMAP for Harmony integration")
  integ_srat_harm <- RunUMAP(integ_srat_harm, reduction = "harmony", dims = 1:ndims)

  # cluster cells
  print("Clustering Harmony integration")
  source("scrnaseq_Leo/utilities/cluster_multiple_res.R")
  integ_srat_harm <- FindNeighbors(integ_srat_harm, reduction = "harmony", dims = 1:ndims) %>%
    cluster_multiple_res(resolutions, "Harmony", out_directory)

  # Normalizing data slot in RNA assay. 
  # Some datasets don't have their data slot in RNA normalized, causing dowstream visualization issues
  integ_srat_harm <- NormalizeData(integ_srat_harm, assay = "RNA")

  if(use_rds){
    # Save integrated dataset as a .rds file
    print("Saving Harmony integration as .rds file")
    file_name = paste(format(Sys.Date(), "%Y%m%d"), "harmony_integ.rds", sep = "_")
    saveRDS(object = integ_srat_harm, file = file.path(out_directory, file_name))
  } else {
    # Save integrated dataset as a .qs file
    print("Saving Harmony integration as .qs file")
    file_name = paste(format(Sys.Date(), "%Y%m%d"), "harmony_integ.qs", sep = "_")
    qsave(integ_srat_harm, file.path(out_directory, file_name))
  }
  
  message("Dataset successfully underwent Harmony integration. The integrated dataset is called: ", file_name)
  
  source("scrnaseq_Leo/pipeline_scripts/visualization.R")
  visualization(integ_srat_harm, ndims = ndims, resolutions = resolutions, integ_type = "Harmony", do_tsne = do_tsne, out_directory = out_directory)

  if(!is.null(log_file)){
    sink(type = "output")
    sink(type = "message")
  }
}