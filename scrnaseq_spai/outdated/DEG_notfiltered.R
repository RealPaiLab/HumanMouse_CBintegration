# Redo of DEG analysis without post filtering.

rm(list=ls())
library(tidyverse)
library(Seurat)

set.seed(1)
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBC1_marker_analysis/Integrated_RL_Hsonly/"
dt <- format(Sys.Date(),"%y%m%d")

outDir <- sprintf("%s/DEG_nofiltered_%s",outDir, dt)
if (!file.exists(outDir)) dir.create(outDir)

cat("Reading dataset\n")
rl_cca <- qs::qread("/home/rstudio/isilon/private/llau/results/integrated/20240618/rl_cca.qs")

ubc_cells <- rl_cca %>%
   subset(subset = species == "human") %>%
   subset(subset = broad_cell_type == "UBC")

cat("Finding intersecting genes")
 ald_mega <- subset(ubc_cells, subset = dataset_name == "Aldinger_full_cerebellum_human")
 sepp_mega <- subset(ubc_cells, subset = dataset_name == "Sepp_full_cerebellum_human")

ald <- ald_mega
sepp <- sepp_mega

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
## # Note that the number of genes in the Aldinger and Sepp datasets are the same as they were merged together before
cat(sprintf("Aldinger: %i total, %i in 1%%\n", 
    length(rownames(ald)),length(ald_genes)))
cat(sprintf("Sepp: %i total, %i in 1%%\n", 
    length(rownames(sepp)),length(sepp_genes)))
cat(sprintf("Intersecting genes = %i",length(intersecting_genes)))
ald <- subset(ald, features = intersecting_genes)
sepp <- subset(sepp, features = intersecting_genes)
cat("Running SCTransform on each set\n")
cat("aldinger\n")
t0 <- Sys.time()
 ald <- SCTransform(ald, vars.to.regress = "CC.Difference", variable.features.n = 5000)
 print(Sys.time()-t0)
cat("Sepp\n")
t0 <- Sys.time()
 sepp <- SCTransform(sepp, vars.to.regress = "CC.Difference", variable.features.n = 5000)
 print(Sys.time()-t0)
cat("Merging UBC\n")
ubc_merged <- merge(x = ald, y = sepp)
cat("PrepSCTFindMarkers\n")
 ubc_merged <- PrepSCTFindMarkers(object = ubc_merged)
cat("Releveling\n")
 # Releveling based on alphabetical then numeric order
 levels_index <- str_order(
   unique(ubc_merged$new_clusters),
   numeric = TRUE
 )
 ubc_merged$new_clusters <- fct_relevel(ubc_merged$new_clusters, levels(ubc_merged$new_clusters)[levels_index])
## # Setting the Idents to be the cluster name because FindConservedMarkers() uses the Ident
Idents(ubc_merged) <- "new_clusters"
# Creating empty list to store markers
cat("Running DEG analysis\n")
ubc_marker_list <- list()
# Find markers for each cluster
for (cluster in levels(Idents(ubc_merged))) {
    print(cluster)
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
cat("merging in DF\n")
# Creating a dataframe from the list of markers
df_ubc <- ubc_marker_list
for (name in names(df_ubc)) {
  # Adding cluster name as a column
  df_ubc[[name]]$cluster <- name
  # Adding gene name as a column
  df_ubc[[name]]$feature <- rownames(df_ubc[[name]])
}
combined_df <- do.call(rbind, df_ubc)
outFile <- sprintf("%s/DEG_%s.txt",outDir,dt)
write.table(combined_df,file=outFile,sep="\t",col=T,row=F,quote=F)

 # }
 #}