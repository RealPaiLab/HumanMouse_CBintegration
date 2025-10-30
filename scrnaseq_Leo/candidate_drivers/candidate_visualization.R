library(argparse)
library(Seurat)
library(qs)
library(tidyverse)
library(patchwork)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # path to list of candidate drivers csv
  "--driver_list",
  default = NULL,
  required = TRUE
)

parser$add_argument(
  # path to out directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--driver_list", "/.mounts/labs/pailab/private/llau/results/integrated/test/drivers.csv",
    "--out_dir", "/.mounts/labs/pailab/private/llau/results/integrated/test"
  ))
} else {
  arg_list <- parser$parse_args()
}

out_dir <- arg_list$out_dir
if (!dir.exists(out_dir )) {
  dir.create(out_dir, recursive = TRUE)
}

drivers <- read.csv(arg_list$driver_list)$drivers

rl_srat <- qread("/.mounts/labs/pailab/private/llau/results/integrated/20240618/rl_cca.qs")

rl_human <- subset(rl_srat, subset = species == "human")
rl_human <- NormalizeData(rl_human, assay = "RNA")
DefaultAssay(rl_human) <- "RNA"
rl_human <- RunUMAP(rl_human, reduction = "pca", dims = 1:25)
human_cells <- length(Cells(rl_human))

tumour_srat <- readRDS("/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds")
tumour_srat <- NormalizeData(tumour_srat, assay = "RNA")
DefaultAssay(tumour_srat) <- "RNA"
tumour_cells <- length(Cells(tumour_srat))

cell_type_plt <- DimPlot(rl_human, reduction = "umap", group.by = "broad_w_ubc_subtypes", label = TRUE, label.size = 4, raster = FALSE, repel = TRUE) +
  ggtitle(paste0("Rhombic Lip Cell Types ", "(n = ", human_cells, ")")) 
filename <- paste0("integ_umap_clusters.png")
#ggsave(filename = filename, plot = cell_type_plt, path = out_dir, 
            #width = 7, height = 7, dpi = 600)

cancer_type_plt <- DimPlot(tumour_srat, reduction = "umap", group.by = "subtype", label = TRUE, label.size = 5, raster = FALSE) +
  ggtitle(paste0("Medulloblastoma Subgroups ", "(n = ", tumour_cells, ")"))
filename <- paste0("integ_umap_clusters.png")
#ggsave(filename = filename, plot = cell_type_plt, path = out_dir, 
            #width = 10, height = 10, dpi = 600)

combined_plot_list <- list()

for (driver in drivers) {
    # Creating out directory
    driver_out <- paste(out_dir, driver, sep = "/")
    if (!dir.exists(driver_out)) {
        dir.create(paste(out_dir, driver, sep = "/"), recursive = TRUE)
    }
    # Feature plot RL
    feature_plot_rl <- FeaturePlot(object = rl_human, features = driver, order = TRUE, min.cutoff = "q10", max.cutoff = "q90", raster = FALSE) +
      ggtitle(paste0(driver, " Rhombic Lip Expression"))
    feature_plot_name <- paste0(driver, "_rl_feature_plot.png")
    #ggsave(filename = feature_plot_name, plot = feature_plot_rl, path = driver_out, 
            #width = 10, height = 10, units = "in", dpi = 600)

    # Feature plot tumour
    feature_plot_tumour <- FeaturePlot(object = tumour_srat, features = driver, order = TRUE, min.cutoff = "q10", max.cutoff = "q90", raster = FALSE) +
      ggtitle(paste0(driver, " Medulloblastoma Expression"))
    feature_plot_name <- paste0(driver, "_tumour_feature_plot.png")
    #ggsave(filename = feature_plot_name, plot = feature_plot_tumour, path = driver_out, 
            #width = 10, height = 10, units = "in", dpi = 600)
    
    dot_plot <- DotPlot(object = rl_human, features = driver, group.by = "broad_w_ubc_subtypes") + theme_classic() + coord_flip() + 
        theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis labels by 45 degrees
        axis.text.y = element_text(size = 12),  # Adjust y-axis label size
        legend.text = element_text(size = 12)) +  # Adjust legend text size
      ggtitle(paste0("Gene Expression of ", driver, " for Rhombic Lip Cell Types")) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
      )
    dot_plot_name <- paste0(driver, "_dot_plot.png")
    #ggsave(filename = dot_plot_name, plot = dot_plot, path = driver_out, 
            #width = 15, height = 5, units = "in", dpi = 600)

    # Fetch expression data and subtype information
    expression_data <- FetchData(tumour_srat, vars = driver)
    subtype_data <- tumour_srat$subtype
    combined_data <- data.frame(Expression = expression_data[[driver]], Subtype = subtype_data)

    # Create the box plot
    expression_plot <- ggplot(combined_data, aes(x = Subtype, y = Expression)) +
        geom_boxplot() +
        labs(title = paste("Expression of", driver, "by Medulloblastoma Subgroup (scRNA-seq)"),
            x = "Subgroups",
            y = "Expression Level") +
        theme_classic()
    box_plot_name <- paste0(driver, "_box_plot.png")
    ggsave(filename = box_plot_name, plot = expression_plot, path = driver_out, 
            width = 10, height = 10, units = "in", dpi = 600)

    combined_plot <- (cell_type_plt + feature_plot_rl)/(dot_plot)/(cancer_type_plt + feature_plot_tumour)
    combined_plot_name <- paste0(driver, "_combined_plot.pdf")
    ggsave(filename = combined_plot_name, plot = combined_plot, path = driver_out, 
            width = 15, height = 15, units = "in", dpi = 300)

    combined_plot_list[[driver]] <- combined_plot
}

pdf(file.path(out_dir, "all_plots.pdf"), width = 15, height = 15)

for (plot in combined_plot_list) {
  print(plot)
}

dev.off()