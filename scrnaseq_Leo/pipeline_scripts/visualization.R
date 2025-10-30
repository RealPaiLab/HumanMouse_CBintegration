#' @import Seurat
#' @import tidyverse

#' Visualization of integrated datasets: UMAP, TSNE, Dot Plot, Feature Plot
#' 
#' @param integ_srat Seurat dataset object 
#' @param ndims Integer indicating number of dimensions to graph
#' @param resolutions Vector of numerical characters that specify the clustering resolutions.
#' @param integ_type String indicating integration type
#' @param do_tsne Boolean indicating to do T-SNE plots
#' @param out_directory String of output directory



visualization <- function(
  integ_srat,
  ndims = 30,
  resolutions,
  integ_type,
  do_tsne = FALSE,
  out_directory
){

  # Determining the number of human and mice cells
  num_cells <- table(integ_srat@meta.data$species)
  human_cell_count <- num_cells["human"]
  mouse_cell_count <- num_cells["mouse"]

  # Making sure this column is a factor
  integ_srat$common_cell_name <- factor(integ_srat$common_cell_name)

  # Plot multiple cluster resolutions
  for (i in resolutions) {
    resolution_name <- paste0("snn_res.", i)
    plt <- DimPlot(integ_srat, reduction = "umap", group.by = resolution_name, label = TRUE, label.size = 3, raster = FALSE) + 
      ggtitle(paste("Resolution:", i))
    filename <- paste0(integ_type, "_", i, "_integ_umap_clusters.png")
    ggsave(filename = filename, plot = plt, path = out_directory, 
            width = 10, height = 10, dpi = 600)
  }


  cols <- brewer.pal(3,"Dark2")
  cols <- alpha(cols, 0.4)
  names(cols) <- levels(integ_srat$species)

  if(do_tsne){
    reduc_types = c("tsne", "umap")
  } else{
    reduc_types = "umap"
  }

  i <- 1
  for (reduc in reduc_types) {
    for (group in c("seurat_clusters", "species", "dataset_name")) {
      # Adding transparency for species plot
      if(group == "species") {
        plt <- DimPlot(integ_srat, reduction = reduc, label = TRUE, repel = TRUE, group.by = group, cols=cols, raster = FALSE)
      } else {
        plt <- DimPlot(integ_srat, reduction = reduc, label = TRUE, repel = TRUE, group.by = group, raster = FALSE)
      }
      # Determining name based column name
      if(group == "seurat_clusters"){
        type = "clusters"
      }
      if(group == "species"){
        type = "species"
      }
      if(group == "dataset_name"){
        type = "dataset"
      }
      name = paste0(integ_type, "_integ_", reduc, "_", group, ".png")
      ggsave(filename = name, plot = plt, path = out_directory, 
            width = 10, height = 10, units = "in", dpi = 600)
      i <- i+1
    }
  }
  

  # Bar plots

  # Species
  source("software/utilities/cluster_barplot.R")
  species_bar_plot <- cluster_barplot(integ_srat, split.by = "species")
  species_bar_plot_name <- paste0(integ_type, "_integ_species_bar_plot.png")
  ggsave(filename = species_bar_plot_name, plot = species_bar_plot, path = out_directory, 
            width = 20, height = 10, units = "in", dpi = 600)
  
  # Common cell name
  common_cell_type_bar_plot <- cluster_barplot(integ_srat, split.by = "common_cell_name")
  common_cell_type_bar_plot_name <- paste0(integ_type, "_integ_common_cell_type_bar_plot.png")
  ggsave(filename = common_cell_type_bar_plot_name, plot = common_cell_type_bar_plot, path = out_directory, 
            width = 20, height = 10, units = "in", dpi = 600)

  # Cell type
  cell_type_bar_plot <- cluster_barplot(integ_srat, split.by = "cell_type")
  cell_type_bar_plot_name <- paste0(integ_type, "_integ_cell_type_bar_plot.png")
  ggsave(filename = cell_type_bar_plot_name, plot = cell_type_bar_plot, path = out_directory, 
            width = 20, height = 10, units = "in", dpi = 600)

  # Dataset
  dataset_bar_plot <- cluster_barplot(integ_srat, split.by = "dataset_name")
  dataset_bar_plot_name <- paste0(integ_type, "_integ_dataset_bar_plot.png")
  ggsave(filename = dataset_bar_plot_name, plot = dataset_bar_plot, path = out_directory, 
            width = 20, height = 10, units = "in", dpi = 600)
  
  

  human_title <- paste0("Human Cell Types (n = ", human_cell_count, ")")
  mouse_title <- paste0("Mouse Cell Types (n = ", mouse_cell_count, ")")

  # Create Dimplot for original cell types
  umap_human_cells <- DimPlot(
    integ_srat,
    cells = rownames(integ_srat@meta.data[integ_srat$species == "human", ]),
    reduction = "umap",
    group.by = "cell_type",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + 
    ggtitle(human_title) + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))

  umap_mouse_cells <- DimPlot(
    integ_srat,
    cells = rownames(integ_srat@meta.data[integ_srat$species == "mouse", ]),
    reduction = "umap",
    group.by = "cell_type",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + 
    ggtitle(mouse_title) + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))

  ggsave(
    filename = paste(integ_type, "_umap_original_celltypes.png", sep = ""),
    plot = umap_human_cells + umap_mouse_cells,
    path = out_directory,
    width = 20,
    height = 20,
    units = "in",
    dpi = 600
  )
  

  cols <- c(brewer.pal(12,"Paired"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
  names(cols) <- levels(integ_srat$common_cell_name)

  # Create Dimplot for common cell types
  umap_human_cells <- DimPlot(
    integ_srat,
    cells = rownames(integ_srat@meta.data[integ_srat$species == "human", ]),
    cols = cols,
    reduction = "umap",
    group.by = "common_cell_name",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + 
    ggtitle(human_title) + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))

  umap_mouse_cells <- DimPlot(
    integ_srat,
    cells = rownames(integ_srat@meta.data[integ_srat$species == "mouse", ]),
    cols = cols,
    reduction = "umap",
    group.by = "common_cell_name",
    label = TRUE,
    raster = FALSE
  ) + 
    ggtitle(mouse_title) + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))

  ggsave(
    filename = paste(integ_type, "_umap_common_celltypes.png", sep = ""),
    plot = umap_human_cells + umap_mouse_cells,
    path = out_directory,
    width = 20,
    height = 10,
    units = "in",
    dpi = 600
  )


  # Genes of interest
  genes <- read.csv("scrnaseq_Leo/utilities/cell_gene_mapping.csv")
  genes <- genes$gene

  DefaultAssay(integ_srat) <- "RNA"
  
  # Dot plot
  dot_plot <- DotPlot(object = integ_srat, features = genes) + theme_classic()
  dot_plot_name <- paste0(integ_type, "_integ_dot_plot.png")
  ggsave(filename = dot_plot_name, plot = dot_plot, path = out_directory, 
            width = 20, height = 10, units = "in", dpi = 600)
  
  # Feature plot
  feature_plot <- FeaturePlot(object = integ_srat, features = genes, order = TRUE, min.cutoff = "q10", max.cutoff = "q90", raster = FALSE)
  feature_plot_name <- paste0(integ_type, "_integ_feature_plot.png")
  ggsave(filename = feature_plot_name, plot = feature_plot, path = out_directory, 
            width = 25, height = 40, units = "in", dpi = 600)

  #Violin plot
  violin_plot <- VlnPlot(object = integ_srat, features = genes, group.by = "seurat_clusters", raster = FALSE)
  violin_plot_name <- paste0(integ_type, "_integ_violin_plot.png")
  ggsave(filename = violin_plot_name, plot = violin_plot, path = out_directory, 
            width = 45, height = 45, units = "in", dpi = 600)
}