## -----------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(qs)
library(RColorBrewer)


## ----eval=FALSE---------------------------------------------------------------
## library(scProportionTest)


## ----eval=FALSE---------------------------------------------------------------
## set.seed(1)


## ----eval=FALSE---------------------------------------------------------------
## setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")
## 
## integ_rl <- qread("/.mounts/labs/pailab/private/llau/results/integrated/20240524/25_pc_without_luo/25_pc_rl.qs")
## 
## # Subsetting out cells
## integ_rl$snn_res.0.4 <- fct_relevel(integ_rl$snn_res.0.4, str_sort(levels(integ_rl$snn_res.0.4), numeric = TRUE))
## integ_ubc <- subset(x = integ_rl, subset = snn_res.0.4 == 2)
## 
## out_directory <- "/.mounts/labs/pailab/private/llau/results/integrated/20240527"
## 
## DefaultAssay(integ_ubc) <- "integrated"


## ----eval=FALSE---------------------------------------------------------------
## # Plots
## resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8)
## # Retrieve UMAP coordinates
## umap_coords <- Embeddings(integ_ubc, "umap")
## 
## # Calculate minimum and maximum values for each dimension
## min_umap <- apply(umap_coords, 2, min) # Min (x,y)
## max_umap <- apply(umap_coords, 2, max) # Max (x,y)
## 
## min_x <- ceiling(min_umap[1])
## min_y <- ceiling(min_umap[2])
## max_x <- ceiling(max_umap[1])
## max_y <- ceiling(max_umap[2])
## 
## buffer_x <- (max_x - min_x)*0.05
## buffer_y <- (max_y - min_y)*0.05
## 
## min_x <- min_x - buffer_x
## min_y <- min_y - buffer_y
## max_x <- max_x + buffer_x
## max_y <- max_y + buffer_y
## 
## source("scrnaseq_Leo/utilities/cluster_multiple_res.R")
## integ_ubc  <- FindNeighbors(integ_ubc, dims = 1:25) %>%
##   cluster_multiple_res(resolutions, 'CCA', out_directory)
## integ_ubc  <- NormalizeData(integ_ubc, assay = "RNA")
## 
## # Plot multiple cluster resolutions
## for (i in resolutions) {
##   resolution_name <- paste0("snn_res.", i)
##   plt <- DimPlot(integ_ubc, reduction = "umap", group.by = resolution_name, label = TRUE, label.size = 3, raster = FALSE, pt.size = 1) +
##     ggtitle(paste("Resolution:", i)) +
##     xlim(min_x, max_x) +       # Set x-axis limits
##     ylim(min_y, max_y)         # Set y-axis limits
##   filename <- paste0(i, "_cluster_umap.png")
##   ggsave(filename = filename, plot = plt, path = out_directory,
##           width = 10, height = 10, dpi = 600)
## }
## 
## qsave(integ_ubc, paste0(out_directory, "/ubc_subset.qs"))


## ----eval=FALSE---------------------------------------------------------------
## # Loading in the UBC cells
## cluster_7 <- read.csv("/.mounts/labs/pailab/private/llau/results/integrated/20240516/fc_subset_analysis/cluster_7_cells.csv")
## 
## cluster_19 <- read.csv("/.mounts/labs/pailab/private/llau/results/integrated/20240516/fc_subset_analysis/cluster_19_cells.csv")
## 
## cluster_20 <- read.csv("/.mounts/labs/pailab/private/llau/results/integrated/20240516/fc_subset_analysis/cluster_20_cells.csv")
## 
## highlight_human_ubc_cells <- list(cluster_19 = cluster_19$x, cluster_20 = cluster_20$x)
## highlight_colours <- c("blue", "green", "red")
## 
## # UBC Cluster 19, 20 DimPlot
## ubc_plt <- DimPlot(integ_ubc, group.by = "species", reduction = "umap", cells.highlight = highlight_human_ubc_cells,
##                   cols.highlight = highlight_colours, label = TRUE, label.size = 3, raster = FALSE) +
##             ggtitle("Cluster 19 and 20 Cells") +
##             xlim(min_x, max_x) +       # Set x-axis limits
##             ylim(min_y, max_y)         # Set y-axis limits
## filename <- "cluster_19_20_cells.png"
## ggsave(filename = filename, plot = ubc_plt, path = out_directory,
##     width = 10, height = 10, dpi = 600)
## 
## highlight_all_human_ubc_cells <- list(cluster_19 = cluster_19$x, cluster_20 = cluster_20$x, cluster_7 = cluster_7$x)
## highlight_colours <- c("red", "blue", "green")
## 
## # UBC Cluster 7, 19, 20 DimPlot
## ubc_plt <- DimPlot(integ_ubc, group.by = "species", reduction = "umap", cells.highlight = highlight_all_human_ubc_cells,
##                   cols.highlight = highlight_colours, label = TRUE, label.size = 3, raster = FALSE) +
##             ggtitle("Cluster 7, 19 and 20 Cells") +
##             xlim(min_x, max_x) +       # Set x-axis limits
##             ylim(min_y, max_y)         # Set y-axis limits
## filename <- "cluster_7_19_20_cells.png"
## ggsave(filename = filename, plot = ubc_plt, path = out_directory,
##         width = 10, height = 10, dpi = 600)
## 
## highlight_colours <- "red"
## 
## # UBC Cluster 7 DimPlot
## ubc_plt <- DimPlot(integ_ubc, group.by = "species", reduction = "umap", cells.highlight = highlight_all_human_ubc_cells$cluster_7,
##                   cols.highlight = highlight_colours, label = TRUE, label.size = 3, raster = FALSE) +
##             ggtitle("Cluster 7") +
##             xlim(min_x, max_x) +       # Set x-axis limits
##             ylim(min_y, max_y)         # Set y-axis limits
## filename <- "cluster_7_cells.png"
## ggsave(filename = filename, plot = ubc_plt, path = out_directory,
##         width = 10, height = 10, dpi = 600)


## ----eval=FALSE---------------------------------------------------------------
## num_cells <- table(integ_ubc@meta.data$species)
## human_cell_count <- num_cells["human"]
## mouse_cell_count <- num_cells["mouse"]
## 
## human_title <- paste0("Human Cell Ages (n = ", human_cell_count, ")")
## mouse_title <- paste0("Mouse Cell Ages (n = ", mouse_cell_count, ")")
## 
## num_bins_human <- 8
## colour_palette_human <- scales::viridis_pal(option = "magma")(num_bins_human)
## 
## num_bins_mouse <- 13
## colour_palette_mouse <- scales::viridis_pal(option = "magma")(num_bins_mouse)
## 
## umap_human_cells <- DimPlot(
##     integ_ubc,
##     cells = rownames(integ_ubc@meta.data[integ_ubc$species == "human", ]),
##     cols = colour_palette_human,
##     reduction = "umap",
##     group.by = "human_age",
##     label = TRUE,
##     repel = TRUE,
##     raster = FALSE,
##     pt.size = 0.5
## ) +
##     ggtitle(human_title) +
##     theme(legend.position = "bottom") +
##     guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2))) +
##     xlim(min_x, max_x) +       # Set x-axis limits
##     ylim(min_y, max_y)         # Set y-axis limits
## 
## umap_mouse_cells <- DimPlot(
##     integ_ubc,
##     cells = rownames(integ_ubc@meta.data[integ_ubc$species == "mouse", ]),
##     cols = colour_palette_mouse,
##     reduction = "umap",
##     group.by = "mouse_age",
##     label = TRUE,
##     repel = TRUE,
##     raster = FALSE,
##     pt.size = 0.5
##   ) +
##     ggtitle(mouse_title) +
##     theme(legend.position = "bottom") +
##     guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2))) +
##     xlim(min_x, max_x) +       # Set x-axis limits
##     ylim(min_y, max_y)         # Set y-axis limits
## 
##   ggsave(
##     filename = "ages_colour_scale_umap.png",
##     plot = umap_human_cells + umap_mouse_cells,
##     path = out_directory,
##     width = 20,
##     height = 15,
##     units = "in",
##     dpi = 600
##   )


## ----eval=FALSE---------------------------------------------------------------
## # 0.1 resolution
## integ_ubc$snn_res.0.1 <- fct_relevel(integ_ubc$snn_res.0.1, str_sort(levels(integ_ubc$snn_res.0.1), numeric = TRUE))
## 
## source("/u/llau/software/mb_scrnaseq/MB_scRNAseq/software/utilities/cluster_barplot.R")
## human_age_bar_plot <- cluster_barplot(integ_ubc, split.by = "human_age", group.by = "snn_res.0.1")
## mouse_age_bar_plot <- cluster_barplot(integ_ubc, split.by = "mouse_age", group.by = "snn_res.0.1")
## bar_plot_name <- "0.1_age_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = human_age_bar_plot + mouse_age_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## dataset_bar_plot <- cluster_barplot(integ_ubc, split.by = "dataset_name", group.by = "snn_res.0.1")
## bar_plot_name <- "0.1_dataset_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = dataset_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## cell_type_bar_plot <- cluster_barplot(integ_ubc, split.by = "cell_type", group.by = "snn_res.0.1")
## bar_plot_name <- "0.1_cell_type_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = cell_type_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## common_cell_type_bar_plot <- cluster_barplot(integ_ubc, split.by = "common_cell_name", group.by = "snn_res.0.1")
## bar_plot_name <- "0.1_common_cell_name_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = common_cell_type_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## species_bar_plot <- cluster_barplot(integ_ubc, split.by = "species", group.by = "snn_res.0.1")
## bar_plot_name <- "0.1_species_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = species_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)


## ----eval=FALSE---------------------------------------------------------------
## cluster_7_metadata <- integ_ubc@meta.data[highlight_all_human_ubc_cells$cluster_7, , drop = FALSE]
## cluster_7_counts <- table(cluster_7_metadata$snn_res.0.1)
## 
## cluster_19_metadata <- integ_ubc@meta.data[highlight_all_human_ubc_cells$cluster_19, , drop = FALSE]
## cluster_19_counts <- table(cluster_19_metadata$snn_res.0.1)
## 
## cluster_20_metadata <- integ_ubc@meta.data[highlight_all_human_ubc_cells$cluster_20, , drop = FALSE]
## cluster_20_counts <- table(cluster_20_metadata$snn_res.0.1)
## 
## all_counts <- data.frame(
##   snn_res_0.1 = rep(names(cluster_7_counts), 3),
##   Count = c(as.vector(cluster_7_counts), as.vector(cluster_19_counts), as.vector(cluster_20_counts)),
##   Cluster = rep(c("Cluster 7", "Cluster 19", "Cluster 20"), each = length(cluster_7_counts))
## )
## 
## cluster_barplot <- ggplot(all_counts, aes(x = snn_res_0.1, y = Count, fill = Cluster)) +
##   geom_bar(stat = "identity", position = "stack", color = "black") +
##   scale_y_continuous(expand = expansion(mult = c(0, .1))) +
##   geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3) +
##   labs(title = "Cluster Counts", x = "snn_res_0.1", y = "Count") +
##   theme_classic()
## 
## bar_plot_name <- "0.1_original_clusters_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = cluster_barplot, path = out_directory,
##        width = 20, height = 10, units = "in", dpi = 600)


## ----eval=FALSE---------------------------------------------------------------
## # 0.1 resolution
## ubc_test_0.1 <- sc_utils(integ_ubc)
## ubc_test_0.1 <- permutation_test(
## 	ubc_test_0.1, cluster_identity = "snn_res.0.1",
## 	sample_1 = "mouse", sample_2 = "human",
## 	sample_identity = "species",
##   n_permutations = 10000
## )
## 
## plt <- permutation_plot(ubc_test_0.1)
## filename <- "0.1_permutation_test.png"
## ggsave(filename = filename, plot = plt, path = out_directory, width = 10, height = 10, dpi = 600)
## ubc_test_0.1@results$permutation

## ----eval=FALSE---------------------------------------------------------------
##    clusters     human     mouse obs_log2FD      pval       FDR boot_mean_log2FD
##      <char>     <num>     <num>      <num>     <num>     <num>            <num>
## 1:        0 0.5555062 0.4817898   0.205399 9.999e-05 9.999e-05        0.2061019
## 2:        1 0.4444938 0.5182102  -0.221374 9.999e-05 9.999e-05       -0.2217122
##    boot_CI_2.5 boot_CI_97.5
##          <num>        <num>
## 1:   0.1436922    0.2684751
## 2:  -0.2836354   -0.1576604


## ----eval=FALSE---------------------------------------------------------------
## # 0.2 resolution
## integ_ubc$snn_res.0.2 <- fct_relevel(integ_ubc$snn_res.0.2, str_sort(levels(integ_ubc$snn_res.0.2), numeric = TRUE))
## 
## source("/u/llau/software/mb_scrnaseq/MB_scRNAseq/software/utilities/cluster_barplot.R")
## human_age_bar_plot <- cluster_barplot(integ_ubc, split.by = "human_age", group.by = "snn_res.0.2")
## mouse_age_bar_plot <- cluster_barplot(integ_ubc, split.by = "mouse_age", group.by = "snn_res.0.2")
## bar_plot_name <- "0.2_age_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = human_age_bar_plot + mouse_age_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## dataset_bar_plot <- cluster_barplot(integ_ubc, split.by = "dataset_name", group.by = "snn_res.0.2")
## bar_plot_name <- "0.2_dataset_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = dataset_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## cell_type_bar_plot <- cluster_barplot(integ_ubc, split.by = "cell_type", group.by = "snn_res.0.2")
## bar_plot_name <- "0.2_cell_type_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = cell_type_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## common_cell_type_bar_plot <- cluster_barplot(integ_ubc, split.by = "common_cell_name", group.by = "snn_res.0.2")
## bar_plot_name <- "0.2_common_cell_name_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = common_cell_type_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## species_bar_plot <- cluster_barplot(integ_ubc, split.by = "species", group.by = "snn_res.0.2")
## bar_plot_name <- "0.2_species_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = species_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)


## ----eval=FALSE---------------------------------------------------------------
## cluster_7_metadata <- integ_ubc@meta.data[highlight_all_human_ubc_cells$cluster_7, , drop = FALSE]
## cluster_7_counts <- table(cluster_7_metadata$snn_res.0.2)
## 
## cluster_19_metadata <- integ_ubc@meta.data[highlight_all_human_ubc_cells$cluster_19, , drop = FALSE]
## cluster_19_counts <- table(cluster_19_metadata$snn_res.0.2)
## 
## cluster_20_metadata <- integ_ubc@meta.data[highlight_all_human_ubc_cells$cluster_20, , drop = FALSE]
## cluster_20_counts <- table(cluster_20_metadata$snn_res.0.2)
## 
## all_counts <- data.frame(
##   snn_res_0.2 = rep(names(cluster_7_counts), 3),
##   Count = c(as.vector(cluster_7_counts), as.vector(cluster_19_counts), as.vector(cluster_20_counts)),
##   Cluster = rep(c("Cluster 7", "Cluster 19", "Cluster 20"), each = length(cluster_7_counts))
## )
## 
## cluster_barplot <- ggplot(all_counts, aes(x = snn_res_0.2, y = Count, fill = Cluster)) +
##   geom_bar(stat = "identity", position = "stack", color = "black") +
##   scale_y_continuous(expand = expansion(mult = c(0, .1))) +
##   geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3) +
##   labs(title = "Cluster Counts", x = "snn_res_0.2", y = "Count") +
##   theme_classic()
## 
## bar_plot_name <- "0.2_original_clusters_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = cluster_barplot, path = out_directory,
##        width = 20, height = 10, units = "in", dpi = 600)


## ----eval=FALSE---------------------------------------------------------------
## # 0.2 resolution
## ubc_test_0.2 <- sc_utils(integ_ubc)
## ubc_test_0.2 <- permutation_test(
## 	ubc_test_0.2, cluster_identity = "snn_res.0.2",
## 	sample_1 = "mouse", sample_2 = "human",
## 	sample_identity = "species",
##   n_permutations = 10000
## )
## 
## plt <- permutation_plot(ubc_test_0.2)
## filename <- "0.2_permutation_test.png"
## ggsave(filename = filename, plot = plt, path = out_directory, width = 10, height = 10, dpi = 600)
## ubc_test_0.2@results$permutation

## ----eval=FALSE---------------------------------------------------------------
##    clusters      human      mouse obs_log2FD      pval       FDR
##      <char>      <num>      <num>      <num>     <num>     <num>
## 1:        0 0.51385801 0.31945890  0.6857395 9.999e-05 9.999e-05
## 2:        1 0.09930339 0.33229275 -1.7425399 9.999e-05 9.999e-05
## 3:        2 0.20157107 0.08671523  1.2169313 9.999e-05 9.999e-05
## 4:        3 0.14287832 0.09746792  0.5517877 9.999e-05 9.999e-05
## 5:        4 0.04238921 0.16406521 -1.9525004 9.999e-05 9.999e-05
##    boot_mean_log2FD boot_CI_2.5 boot_CI_97.5
##               <num>       <num>        <num>
## 1:        0.6863067   0.6049694    0.7704627
## 2:       -1.7427373  -1.8712664   -1.6166691
## 3:        1.2185864   1.0410970    1.4113252
## 4:        0.5533755   0.3779327    0.7400280
## 5:       -1.9531279  -2.1590142   -1.7557423


## ----eval=FALSE---------------------------------------------------------------
## # 0.3 resolution
## integ_ubc$snn_res.0.2 <- fct_relevel(integ_ubc$snn_res.0.3, str_sort(levels(integ_ubc$snn_res.0.3), numeric = TRUE))
## 
## source("/u/llau/software/mb_scrnaseq/MB_scRNAseq/software/utilities/cluster_barplot.R")
## human_age_bar_plot <- cluster_barplot(integ_ubc, split.by = "human_age", group.by = "snn_res.0.3")
## mouse_age_bar_plot <- cluster_barplot(integ_ubc, split.by = "mouse_age", group.by = "snn_res.0.3")
## bar_plot_name <- "0.3_age_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = human_age_bar_plot + mouse_age_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## dataset_bar_plot <- cluster_barplot(integ_ubc, split.by = "dataset_name", group.by = "snn_res.0.3")
## bar_plot_name <- "0.3_dataset_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = dataset_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## cell_type_bar_plot <- cluster_barplot(integ_ubc, split.by = "cell_type", group.by = "snn_res.0.3")
## bar_plot_name <- "0.3_cell_type_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = cell_type_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## common_cell_type_bar_plot <- cluster_barplot(integ_ubc, split.by = "common_cell_name", group.by = "snn_res.0.3")
## bar_plot_name <- "0.3_common_cell_name_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = common_cell_type_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)
## 
## species_bar_plot <- cluster_barplot(integ_ubc, split.by = "species", group.by = "snn_res.0.3")
## bar_plot_name <- "0.3_species_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = species_bar_plot, path = out_directory,
##           width = 20, height = 10, units = "in", dpi = 600)


## ----eval=FALSE---------------------------------------------------------------
## cluster_7_metadata <- integ_ubc@meta.data[highlight_all_human_ubc_cells$cluster_7, , drop = FALSE]
## cluster_7_counts <- table(cluster_7_metadata$snn_res.0.3)
## 
## cluster_19_metadata <- integ_ubc@meta.data[highlight_all_human_ubc_cells$cluster_19, , drop = FALSE]
## cluster_19_counts <- table(cluster_19_metadata$snn_res.0.3)
## 
## cluster_20_metadata <- integ_ubc@meta.data[highlight_all_human_ubc_cells$cluster_20, , drop = FALSE]
## cluster_20_counts <- table(cluster_20_metadata$snn_res.0.3)
## 
## all_counts <- data.frame(
##   snn_res_0.3 = rep(names(cluster_7_counts), 3),
##   Count = c(as.vector(cluster_7_counts), as.vector(cluster_19_counts), as.vector(cluster_20_counts)),
##   Cluster = rep(c("Cluster 7", "Cluster 19", "Cluster 20"), each = length(cluster_7_counts))
## )
## 
## cluster_barplot <- ggplot(all_counts, aes(x = snn_res_0.3, y = Count, fill = Cluster)) +
##   geom_bar(stat = "identity", position = "stack", color = "black") +
##   scale_y_continuous(expand = expansion(mult = c(0, .1))) +
##   geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3) +
##   labs(title = "Cluster Counts", x = "snn_res_0.3", y = "Count") +
##   theme_classic()
## 
## bar_plot_name <- "0.3_original_clusters_bar_plot.png"
## ggsave(filename = bar_plot_name, plot = cluster_barplot, path = out_directory,
##        width = 20, height = 10, units = "in", dpi = 600)


## ----eval=FALSE---------------------------------------------------------------
## # 0.3 resolution
## ubc_test_0.3 <- sc_utils(integ_ubc)
## ubc_test_0.3 <- permutation_test(
## 	ubc_test_0.3, cluster_identity = "snn_res.0.3",
## 	sample_1 = "mouse", sample_2 = "human",
## 	sample_identity = "species",
##   n_permutations = 10000
## )
## 
## plt <- permutation_plot(ubc_test_0.3)
## filename <- "0.3_permutation_test.png"
## ggsave(filename = filename, plot = plt, path = out_directory, width = 10, height = 10, dpi = 600)
## 
## write.csv(
##   x = ubc_test_0.3@results$permutation,
##   file = file.path(out_directory, "cca_ubc_proportion_test.csv"),
##   row.names = FALSE
## )

## ----eval=FALSE---------------------------------------------------------------
##    clusters      human      mouse obs_log2FD       pval         FDR
##      <char>      <num>      <num>      <num>      <num>       <num>
## 1:        0 0.24722099 0.21262574  0.2174850 0.00019998 0.000199980
## 2:        1 0.26382096 0.10613944  1.3135984 0.00009999 0.000119988
## 3:        2 0.21253891 0.08810267  1.2704693 0.00009999 0.000119988
## 4:        3 0.09915518 0.33055845 -1.7371453 0.00009999 0.000119988
## 5:        4 0.14969616 0.10926119  0.4542562 0.00009999 0.000119988
## 6:        5 0.02756781 0.15331252 -2.4754191 0.00009999 0.000119988
##    boot_mean_log2FD boot_CI_2.5 boot_CI_97.5
##               <num>       <num>        <num>
## 1:        0.2184699   0.1031931    0.3331119
## 2:        1.3151291   1.1500711    1.4841063
## 3:        1.2720907   1.0952528    1.4563973
## 4:       -1.7378439  -1.8667054   -1.6108455
## 5:        0.4566767   0.2872559    0.6339961
## 6:       -2.4816889  -2.7286280   -2.2474339


## ----eval=FALSE---------------------------------------------------------------
## library(tidyverse)
## library(Seurat)
## library(qs)
## library(clustree)

## ----eval=FALSE---------------------------------------------------------------
## integ_ubc <- qread("/.mounts/labs/pailab/private/llau/results/integrated/20240527/ubc_subset.qs")
## out_directory <- "/.mounts/labs/pailab/private/llau/results/integrated/20240527"
## tree <- clustree(integ_ubc, prefix = "snn_res.")
## tree_plot_name <- "tree_plot.png"
## ggsave(filename = tree_plot_name, plot = tree, path = out_directory,
##         width = 15, height = 25, units = "in", dpi = 600)


## ----eval=FALSE---------------------------------------------------------------
## print(sessionInfo())
## 
## R version 4.2.2 (2022-10-31)
## Platform: x86_64-conda-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.6 LTS
## 
## Matrix products: default
## BLAS/LAPACK: /u/llau/.conda/envs/scproportiontest/lib/libopenblasp-r0.3.27.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base
## 
## other attached packages:
##  [1] RColorBrewer_1.1-3          qs_0.25.7
##  [3] SeuratObject_5.0.2          Seurat_4.3.0.1
##  [5] forcats_1.0.0               stringr_1.5.1
##  [7] dplyr_1.1.4                 purrr_1.0.2
##  [9] readr_2.1.5                 tidyr_1.3.1
## [11] tibble_3.2.1                ggplot2_3.5.1
## [13] tidyverse_1.3.2             scProportionTest_0.0.0.9000
## 
## loaded via a namespace (and not attached):
##   [1] googledrive_2.1.1      Rtsne_0.17             colorspace_2.1-0
##   [4] deldir_2.0-4           ggridges_0.5.6         fs_1.6.4
##   [7] spatstat.data_3.0-4    farver_2.1.2           leiden_0.4.3.1
##  [10] listenv_0.9.1          ggrepel_0.9.5          lubridate_1.9.3
##  [13] fansi_1.0.6            xml2_1.3.3             codetools_0.2-20
##  [16] splines_4.2.2          polyclip_1.10-6        spam_2.10-0
##  [19] jsonlite_1.8.8         broom_1.0.6            ica_1.0-3
##  [22] cluster_2.1.6          dbplyr_2.5.0           png_0.1-8
##  [25] uwot_0.1.16            shiny_1.8.1.1          sctransform_0.4.1
##  [28] spatstat.sparse_3.0-3  compiler_4.2.2         httr_1.4.7
##  [31] backports_1.4.1        Matrix_1.6-5           fastmap_1.2.0
##  [34] lazyeval_0.2.2         gargle_1.5.2           cli_3.6.2
##  [37] later_1.3.2            htmltools_0.5.8.1      tools_4.2.2
##  [40] igraph_1.4.2           dotCall64_1.1-1        gtable_0.3.5
##  [43] glue_1.7.0             RANN_2.6.1             reshape2_1.4.4
##  [46] Rcpp_1.0.12            scattermore_1.2        cellranger_1.1.0
##  [49] vctrs_0.6.5            spatstat.explore_3.2-6 nlme_3.1-164
##  [52] progressr_0.14.0       lmtest_0.9-40          spatstat.random_3.2-3
##  [55] globals_0.16.3         rvest_1.0.4            timechange_0.3.0
##  [58] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.4
##  [61] irlba_2.3.5.1          googlesheets4_1.1.1    goftest_1.2-3
##  [64] future_1.33.2          MASS_7.3-60.0.1        zoo_1.8-12
##  [67] scales_1.3.0           ragg_1.2.5             hms_1.1.3
##  [70] promises_1.3.0         spatstat.utils_3.0-4   parallel_4.2.2
##  [73] reticulate_1.37.0      pbapply_1.7-2          gridExtra_2.3
##  [76] stringi_1.7.12         systemfonts_1.1.0      rlang_1.1.3
##  [79] pkgconfig_2.0.3        matrixStats_1.3.0      lattice_0.22-6
##  [82] ROCR_1.0-11            tensor_1.5             labeling_0.4.3
##  [85] patchwork_1.2.0        htmlwidgets_1.6.4      cowplot_1.1.3
##  [88] tidyselect_1.2.1       parallelly_1.37.1      RcppAnnoy_0.0.22
##  [91] plyr_1.8.9             magrittr_2.0.3         R6_2.5.1
##  [94] generics_0.1.3         DBI_1.2.2              withr_3.0.0
##  [97] haven_2.5.4            pillar_1.9.0           fitdistrplus_1.1-11
## [100] survival_3.6-4         abind_1.4-5            sp_2.1-4
## [103] future.apply_1.11.2    crayon_1.5.2           modelr_0.1.11
## [106] KernSmooth_2.23-24     utf8_1.2.4             RApiSerialize_0.1.2
## [109] spatstat.geom_3.2-9    plotly_4.10.4          tzdb_0.4.0
## [112] readxl_1.4.3           grid_4.2.2             data.table_1.15.2
## [115] reprex_2.1.0           digest_0.6.35          xtable_1.8-4
## [118] httpuv_1.6.15          textshaping_0.3.6      RcppParallel_5.1.6
## [121] stringfish_0.15.7      munsell_0.5.1          viridisLite_0.4.2
## 


## ----eval=FALSE---------------------------------------------------------------
## write.csv(srat[[]], file = ...)

