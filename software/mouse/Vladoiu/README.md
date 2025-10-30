# List of scripts

These scripts were used to analyze the Vladoiu single-cell data.

## `isolate_RL.R`

Script to subset the RL cells from the entire cerebellum.

## `load_vladoiu.R`

Functions to create a Seurat object from the raw count matrices.

## `make_integrated_seurat.R`

Script to integrate all the counts from the different timepoints (E10 to P14) using canonical correlation analysis. Not used for downstream analysis since this resulted in cells from different timepoints clustering together.

## `make_merged_seurat.R`

Script to merge all the counts from the different timepoints (E10 to P14).

## `plot_cells_by_timepoint.R`

Script to plot the number of cells from each timepoint.

## `plot_gene_expression.R`

Creates and saves plots, such as `FeaturePlot` and `VlnPlot`, that show the expression of selected genes.

## `visualize_clusters.R`

Script to perform dimensionality reduction and clustering on the full Seurat object followed by visualization of the results.
