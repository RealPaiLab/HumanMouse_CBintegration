 
 library(tidyverse)
 library(Seurat)
 library(qs)
 library(RColorBrewer)
 library(patchwork)
 library(ComplexHeatmap)
 library(circlize)
 library(ggrepel)

 set.seed(1)
 out_directory <- "/.mounts/labs/pailab/private/llau/results/integrated/20240715"
 if (!dir.exists(out_directory)) {
   dir.create(out_directory, recursive = TRUE)
 }
 rl_cca <- qread("/.mounts/labs/pailab/private/llau/results/integrated/20240618/rl_cca.qs")
 
 ubc_cells <- rl_cca %>%
   subset(subset = species == "human") %>%
   subset(subset = broad_cell_type == "UBC")
 
 ald <- subset(ubc_cells, subset = dataset_name == "Aldinger_full_cerebellum_human")
 sepp <- subset(ubc_cells, subset = dataset_name == "Sepp_full_cerebellum_human")
 
 # Get counts slot from Seurat object
 counts_ald <- GetAssayData(object = ald, slot = "counts")
 counts_sepp <- GetAssayData(object = sepp, slot = "counts")
 
 # If a count for a gene in a cell is greater than 0, set as TRUE (= 1)
 nonzero_ald <- counts_ald > 0
 nonzero_sepp <- counts_sepp > 0
 
 # If 1% or more cells are TRUE, keep the gene. Each TRUE value = 1. Taking the sum of all the cells for that gene
 keep_ald_genes <- Matrix::rowSums(nonzero_ald) >= (0.01*length(Cells(ald)))
 keep_sepp_genes <- Matrix::rowSums(nonzero_sepp) >= (0.01*length(Cells(sepp)))
 
 # Get the genes names that we are keeping
 ald_genes <- rownames(ald)[keep_ald_genes]
 sepp_genes <- rownames(sepp)[keep_sepp_genes]
 
 intersecting_genes <- intersect(ald_genes, sepp_genes)
 
 # Note that the number of genes in the Aldinger and Sepp datasets are the same as they were merged together before
 > print(length(rownames(ald)))
 [1] 52672
 
 > print(length(rownames(sepp)))
 [1] 52672
 
 > print(length(ald_genes))
 [1] 9882
 
 > print(length(sepp_genes))
 [1] 12658
 
 > print(length(intersecting_genes))
 [1] 8122
 
 ald <- subset(ald, features = intersecting_genes)
 sepp <- subset(sepp, features = intersecting_genes)
 
 ald <- SCTransform(ald, vars.to.regress = "CC.Difference", variable.features.n = 5000)
 sepp <- SCTransform(sepp, vars.to.regress = "CC.Difference", variable.features.n = 5000)
 
 ubc_merged <- merge(x = ald, y = sepp)
 
 ubc_merged <- PrepSCTFindMarkers(object = ubc_merged)
 
 # Releveling based on alphabetical then numeric order
 levels_index <- str_order(
   unique(ubc_merged$new_clusters),
   numeric = TRUE
 )
 ubc_merged$new_clusters <- fct_relevel(ubc_merged$new_clusters, levels(ubc_merged$new_clusters)[levels_index])
 
 # Setting the Idents to be the cluster name because FindConservedMarkers() uses the Ident
 Idents(ubc_merged) <- "new_clusters"
 
 # Creating empty list to store markers
 ubc_marker_list <- list()
 
 # Find markers for each cluster
 for (cluster in levels(Idents(ubc_merged))) {
   ubc_markers <- FindConservedMarkers(
     object = ubc_merged,
     ident.1 = cluster,
     grouping.var = "dataset_name",
     logfc.threshold = 0,
     min.pct = 0.01,
     only.pos = FALSE,
     assay = "SCT"
   )
   # Append markers to list
   ubc_marker_list[[cluster]] <- ubc_markers
 }
 
 # Creating a dataframe from the list of markers
 df_ubc <- ubc_marker_list
 for (name in names(df_ubc)) {
   # Adding cluster name as a column
   df_ubc[[name]]$cluster <- name
   # Adding gene name as a column
   df_ubc[[name]]$feature <- rownames(df_ubc[[name]])
 }
 
 # Stack into one df
 combined_df <- do.call(rbind,df_ubc) # %>%