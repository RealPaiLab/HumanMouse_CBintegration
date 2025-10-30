# Redo of DEG analysis without post filtering.

rm(list=ls())
library(tidyverse)
library(Seurat)

set.seed(1)
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/"
dt <- format(Sys.Date(),"%y%m%d")

outDir <- sprintf("%s/DEG_nofiltered_%s",outDir, dt)
if (!file.exists(outDir)) dir.create(outDir)

inFile <- sprintf("%s/RLlineage_only/241107/UBC_withClusterAssignments_241107.qs",outDir)

cat("Reading dataset\n")
srat <- qs::qread(inFile)

mdata <- srat@meta.data

clusterColumn <- "integrated_snn_res.0.2"
Idents(srat) <- clusterColumn

stop("filter for genes that are expressed in 1% of Aldinger and Sepp.")
stop("do we want to filter for low-expressing genes?")

srat <- PrepSCTFindMarkers(srat)

ubc_marker_list<- list()
for (cluster in levels(Idents(srat))) {
    print(cluster)
    ubc_markers <- FindMarkers(
        object = srat,
        ident.1 = cluster,
        logfc.threshold = 0,
        min.pct = 0.01,
        assay = "SCT",
        recorrect_umi = FALSE
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
    outFile <- sprintf("%s/UBC_DEG_%s_%s.txt",outDir,clusterColumn, dt)
    write.table(combined_df,file=outFile,sep="\t",col=T,row=F,quote=F)

 # }
 #}