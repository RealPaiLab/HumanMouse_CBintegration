#' Take unprocessed Seurat object containing only the counts data and processes
#' it, including normalization, dimensional reduction, and visualization.
#' 
#' @param object Unprocessed Seurat object
#' 
#' @returns Processed Seurat object
#' 
process_seurat <- function(
  object,
  out_path,
  ndims = 30,
  SCTransform_args = list(vars.to.regress = "CC.Difference"),
  DimHeatmap_args = list(),
  FindClusters_args = list(),
  RunTSNE_args = list(),
  RunUMAP_args = list()
) {
  # set up folder directories
  pca_dir <- file.path(out_path, "pca")
  dir.create(pca_dir)
  tsne_dir <- file.path(out_path, "tsne")
  dir.create(tsne_dir)
  umap_dir <- file.path(out_path, "umap")
  dir.create(umap_dir)
  
  # run sctransform to normalize and find variable features
  object <- do.call(
    SCTransform,
    args = c(list(object = object), SCTransform_args)
  )
  
  # get variable features
  # var_feat <- VariableFeatures(object)
  
  # run PCA
  object <- RunPCA(object)
  
  # save heatmap
  heatmap <- do.call(
    DimHeatmap,
    args = c(list(object = object), DimHeatmap_args)
  )
  ggsave(
    filename = "heatmap.pdf",
    plot = heatmap,
    path = pca_dir
  )
  
  # save elbow plot
  elbow_plot <- ElbowPlot(object, ndims = 50)
  ggsave(
    filename = "elbow_plot.pdf",
    plot = elbow_plot,
    path = pca_dir,
    width = 6,
    height = 4,
    units = "in"
  )
  
  # cluster cells
  object <- FindNeighbors(object, dims = 1:ndims)
  object <- do.call(
    FindClusters,
    args = c(list(object = object), FindClusters_args)
  )
  
  # run t-SNE
  object <- do.call(
    RunTSNE,
    args = c(list(object = object, dims = 1:ndims), RunTSNE_args)
  )
  tsne_plot <- DimPlot(object, reduction = "tsne", label = TRUE, repel = TRUE) + 
    NoLegend()
  ggsave(
    "clusters.pdf",
    plot = tsne_plot,
    path = tsne_dir,
    width = 5,
    height = 5,
    units = "in"
  )
  
  # run UMAP
  object <- do.call(
    RunUMAP,
    args = c(list(object = object, dims = 1:ndims), RunUMAP_args)
  )
  umap_plot <- DimPlot(object, reduction = "umap", label = TRUE, repel = TRUE) + 
    NoLegend()
  ggsave(
    "clusters.pdf",
    plot = umap_plot,
    path = umap_dir,
    width = 5,
    height = 5,
    units = "in"
  )
  
  return(object)
}


#' Takes list of Seurat objects and integrates them. Returns an integrated
#' Seurat object.
#' 
integrate_seurat <- function(
  object.list,
  nfeatures = 2000,
  vars.to.regress = NULL
) {
  # normalize and find variable features
  object.list <- lapply(
    X = object.list,
    FUN = function(X, nfeatures) {
      X <- NormalizeData(X) %>% 
        FindVariableFeatures(., selection.method = "vst", nfeatures = nfeatures)
    },
    nfeatures = nfeatures
  )
  
  # perform integration
  features <- SelectIntegrationFeatures(object.list, nfeatures = nfeatures)
  anchors <- FindIntegrationAnchors(object.list, anchor.features = features)
  object.integrated <- IntegrateData(anchors)
  object.integrated <- ScaleData(object.integrated, vars.to.regress = vars.to.regress)
  
  return(object.integrated)
}
