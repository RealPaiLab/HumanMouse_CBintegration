library(tidyverse)
library(Seurat)
library(qs)
library(RColorBrewer)
library(patchwork)
library(ggrepel)
library(ComplexHeatmap)

set.seed(1)
out_directory <- "/.mounts/labs/pailab/private/llau/results/human/20241108/rl"
if (!dir.exists(out_directory)) {
  dir.create(out_directory, recursive = TRUE)
}
dataset_srat <- readRDS("/.mounts/labs/pailab/private/llau/data/Aldinger_2021/glutamatergic_dev_Liam.RDS")

dataset_srat <- SCTransform(dataset_srat, vars.to.regress = c("CC.Difference", "sex"), variable.features.n = 5000) %>%
  RunPCA()

# Is ScaleData() needed for the heatmaps?
dataset_srat <- ScaleData(dataset_srat,assay="SCT",features=rownames(dataset_srat))

plt <- ElbowPlot(dataset_srat, ndims = 50, reduction = "pca") + theme_classic()
ggsave(filename = "dataset_srat_elbow_plot.png", 
  plot = plt, path = out_directory, width = 7, 
  height = 7, units = "in", dpi = 600)

dataset_srat <- RunUMAP(dataset_srat, dims = 1:25) %>% 
  FindNeighbors(dims = 1:25)

# Removing previously generated clusters
pattern <- paste0(".*_snn_res.")
matching_column <- grep(pattern, colnames(dataset_srat@meta.data), value = TRUE)
for(column in matching_column){
    dataset_srat[[column]] <- NULL
}

source("/u/llau/software/mb_scrnaseq/MB_scRNAseq/scrnaseq_Leo/utilities/cluster_multiple_res.R")

# The last two parameters for this function is just denotating the integration type and the out directory which we don't need to use here
resolutions <- seq(0.2, 0.8, by = 0.2)
dataset_srat <- cluster_multiple_res(
  dataset_srat,
  resolutions,
  "",
  ""
)

plt <- DimPlot(dataset_srat, reduction = "umap", label = TRUE, 
repel = TRUE, group.by = "new_cell_type")
ggsave(filename = "ald_original_cell_annotation_umap.png", 
  plot = plt, path = out_directory, width = 7, 
  height = 7, units = "in", dpi = 600)

for (i in resolutions) {
  out_folder <- paste0(out_directory, "/", i, "_rl_res")
  print(out_folder)
  pattern <- paste0("snn_res.", i)
  if (!dir.exists(out_folder)) {
    dir.create(out_folder, recursive = TRUE)
  }
  dt <- format(Sys.Date(),"%y%m%d")

  #Plotting the UMAP of all clusters in resolution
  plt <- DimPlot(dataset_srat, reduction = "umap", group.by = pattern, label = TRUE, label.size = 3) + ggtitle(paste("Resolution:", i))

  ggsave(filename = paste0("dataset_srat_", i, "_clustering_umap.png"), 
  plot = plt, path = out_folder, width = 7, 
  height = 7, units = "in", dpi = 600)


  # Plotting individual clusters
  plist <- list()

    for (k in levels(dataset_srat@meta.data[[pattern]])) {
      dataset_srat$cur_cluster <- dataset_srat@meta.data[[pattern]] == k
      plist[[k]] <- DimPlot(dataset_srat  , reduction="umap", 
          cols=c("grey90","#aa0000"),
          group.by="cur_cluster",
          order=TRUE,
          label=FALSE) + ggtitle(k) + 
          xlab("UMAP 1") + ylab("UMAP 2") + 
          theme(plot.title=element_text(size=9)) + NoLegend()
    }
    p <- wrap_plots(plist, ncol = 3)
    ggsave(filename = paste0("dataset_srat_", i, "_individual_clusters_umap.png"), 
      plot = p, path = out_folder, width = 9, 
      height = ceiling(length(levels(dataset_srat@meta.data[[pattern]]))/3)*3, units = "in", dpi = 600)
  
  # Plotting dot plot
  genes <- c("MKI67","WLS", "EOMES", "LMX1A", "OTX2", "RBFOX3", "ATOH1","PAX6", "RELN")
  Idents(dataset_srat) <- pattern

  # Note: Plotting the SCT assay
  plt <- DotPlot(object = dataset_srat, assay = "SCT", features = genes) + 
    theme_classic() 
  ggsave(filename = paste0("dataset_srat_", i, "_dot_plot.png"), 
    plot = plt, path = out_folder, width = 10, 
    height = 5, units = "in", dpi = 600)

  # Differential expression

  # Get counts slot from Seurat object
  counts <- GetAssayData(object = dataset_srat, slot = "counts")

  # If a count for a gene in a cell is greater than 0, set as TRUE (= 1)
  nonzero_ald <- counts > 0
  # If 1% or more cells are TRUE, keep the gene. Each TRUE value = 1. Taking the sum of all the cells for that gene
  keep_ald_genes <- Matrix::rowSums(nonzero_ald) >= (0.01*length(Cells(dataset_srat)))
  
  # Get the genes names that we are keeping
  ald_genes <- rownames(dataset_srat)[keep_ald_genes]

  ## # Note that the number of genes in the Aldinger and Sepp datasets are the same as they were merged together before
  cat(sprintf("Aldinger: %i total, %i in 1%%\n", 
      length(rownames(dataset_srat)),length(ald_genes)))

  dataset_srat <- subset(dataset_srat, features = ald_genes)
 
  cat("PrepSCTFindMarkers\n")
  dataset_srat <- PrepSCTFindMarkers(object = dataset_srat)
  cat("Releveling\n")
  # Releveling based on alphabetical then numeric order
  levels_index <- str_order(
    unique(dataset_srat@meta.data[[pattern]]),
    numeric = TRUE
  )
  dataset_srat@meta.data[[pattern]] <- fct_relevel(dataset_srat@meta.data[[pattern]], levels(dataset_srat@meta.data[[pattern]][levels_index]))

  # Setting the Idents to be the cluster name because FindMarkers() uses the Ident
  Idents(dataset_srat) <- pattern

  cat("Running DEG analysis\n")
  marker_list <- list()
  # Find markers for each cluster
  for (cluster in levels(Idents(dataset_srat))) {
      print(cluster)
      markers <- FindMarkers(
        object = dataset_srat,
        ident.1 = cluster,
        logfc.threshold = 0,
        min.pct = 0.01,
        only.pos = FALSE,
        assay = "SCT"
      )
    # Append markers to list
    marker_list[[cluster]] <- markers
  }
  cat("merging in DF\n")
  # Creating a dataframe from the list of markers
  df_cells <- marker_list
  for (name in names(df_cells)) {
    # Adding cluster name as a column
    df_cells[[name]]$cluster <- name
    # Adding gene name as a column
    df_cells[[name]]$feature <- rownames(df_cells[[name]])
  }
  combined_df <- do.call(rbind, df_cells)
  write.csv(combined_df, file = paste0(out_folder, "/", i, "_markers.csv"), row.names = FALSE)

    #Plotting the heatmaps

    DEGresults <- combined_df 

    # create top annotation rows
  dset <- dataset_srat@meta.data[,pattern]
  dset <- as.matrix(table(dset))
  dset <- t(apply(dset,1,function(x) x/sum(x)))
  ages <- as.matrix(table(dataset_srat@meta.data[,c(pattern,"age")]))
  ages <- t(apply(ages,1,function(x) x/sum(x)))
  sex <- as.matrix(table(dataset_srat@meta.data[,c(pattern,"sex")]))
  sex <- t(apply(sex,1,function(x) x/sum(x)))
  pal1 <- RColorBrewer::brewer.pal(ncol(ages),"Spectral") # age
  pal2 <- c("pink","blue") # sex

  # reorder all
  ages <- ages[levels_index,]
  sex <- sex[levels_index,]

  topannot <- HeatmapAnnotation(
      age = anno_barplot(ages, 
          gp = gpar(fill = pal1, border="white",lty="solid"),
          bar_width = 1, 
          height = unit(6, "cm"),
          axis=FALSE
      ),
      sex = anno_barplot(sex,
          gp = gpar(fill = pal2,border="white"),
          bar_width = 1,
          height = unit(2,"cm"),
          axis=FALSE
      )
  )
  lgd_list <- list(
      Legend(labels=colnames(ages),title="age",legend_gp=gpar(fill=pal1)),
      Legend(labels=colnames(sex),title="sex",legend_gp=gpar(fill=pal2))
  )

  orig_dat <- read.csv(DEGresults, header = TRUE, stringsAsFactors = FALSE)

  # -----------------------------------------------------------------
  # Generate heatmaps for different cutoffs

  cutoffs <- c(1,1.5,2)
  for (cutoff in cutoffs) {
      message("---------------------------------")
      message(sprintf("cutoff = %1.2f", cutoff))

      topGenes <- list()
      topData <- list()
      message("collecting top genes per RL cluster")
      for  (cl in levels_index) {
          dat <- subset(orig_dat,cluster %in% cl)
          dat <- subset(dat, 
              p_val < 0.05 
          )
          # Taking p_val, avg_log2FC, and feature columns
          dat <- dat[,c(1,2,7)]
          dat <- dat[order(dat$avg_log2FC,decreasing=TRUE),]
        # message(cl)
          #print(head(dat[,c("avgLogFC","feature")]))
          dat <- subset(dat,
              avg_log2FC > cutoff)
          nr <- min(100,nrow(dat))
          if (nr > 0) {
              topGenes[[cl]] <- dat$feature[1:nr]
              topData[[cl]] <- cbind(rep(cutoff,nr),rep(cl,nr),dat[1:nr,])
          }
          
      }
      topData <- do.call("rbind",topData)
      colnames(topData)[1:2] <- c("avgLogFC_cutoff","topGene_cluster")
      write.table(topData,
          file=sprintf("%s/RL.topGenes.cutoff%1.2f.%s.txt",
              out_folder,cutoff,dt),
          sep="\t",col=TRUE,row=FALSE,quote=FALSE
      )

      message("\tplotting heatmap at cellular resolution")
      feat <- c(unlist(topGenes),
              "SOX4","SOX11","FOXP2","EYS","CRX","NRL")
      
      hm <- DoHeatmap(
          dataset_srat,
          features = feat,
          group.by = pattern,
          size=3
      ) + theme(axis.text.y=element_blank())

      hm <- hm + ggtitle(sprintf("RL markers (Log2FC > %1.2f)",cutoff))
      ggsave(hm, file=sprintf("%s/heatmap_RLmarkers_cutoff%1.2f.%s.png",
          out_folder,cutoff, dt)
      )

      message("plot heatmap at cluster resolution")
      xpr <- as.matrix(AverageExpression(dataset_srat,
          features=feat,
          assay="SCT",
          group.by=pattern
      )[[1]])
      colnames(xpr) <- sub("-","_",colnames(xpr))
      xpr <- xpr[,levels_index]

      quants <- quantile(xpr,c(0.1,0.95))
      col_fun <- circlize::colorRamp2(quants,c("#000000","#FFFF00"))

      annoGenes <- c("CRX","SOX4","SOX11","EYS","NRL")
      message("Genes to annotate")
      for (nm in names(topGenes)) {
          cur <- topGenes[[nm]]; ln <- length(cur)
          cur <- cur[1:min(length(cur),5)]
          annoGenes <- c(annoGenes, cur)
          message(sprintf("%s: %i in set; selected = { %s }", 
              nm, ln, paste(cur,collapse=",")))
      }

      # add top genes per cluster on the right
      ha <- rowAnnotation(foo = anno_mark(at = match(annoGenes,rownames(xpr)),
          labels = annoGenes))
      hm2 <- Heatmap(xpr, name="AverageExpression",
              cluster_rows=FALSE,
              cluster_columns=FALSE,
              col = col_fun,
              show_row_names = FALSE,
              column_title=sprintf("RL markers (log2FC > %1.2f), unscaled",
                  cutoff), 
          top_annotation = topannot,
          right_annotation = ha
      )

      png(sprintf("%s/heatmap_AveRLmarkers_cutoff%1.2f.unscaled.%s.png",
          out_folder, cutoff,dt),
          height=11,width=7,unit="in",res=72
      )
      draw(hm2, annotation_legend_list = lgd_list)
      dev.off()

      message("now plot scaled")
      unscaled <- xpr
      xpr <- t(scale(t(xpr)))
      col_fun <- circlize::colorRamp2(
          c(-2,0,2),c("#FF00FF","#000000", "#FFFF00")
        #c(0, max(xpr)), c("#000000","#FFFF00")
      )

      # add top genes per cluster on the right
      ha <- rowAnnotation(foo = anno_mark(at = match(annoGenes,rownames(xpr)),
          labels = annoGenes))
      hm3 <- Heatmap(xpr, name="AverageExpression",
              cluster_rows=FALSE,
              cluster_columns=FALSE,
              col = col_fun,
              column_title=sprintf("RL markers (log2FC > %1.2f), scaled",
                  cutoff),
              show_row_names = FALSE,
              top_annotation = topannot,
              right_annotation = ha
      )

      png(sprintf("%s/heatmap_AveRLmarkers_cutoff%1.2f.scaled.%s.png",
          out_folder, cutoff,dt),
          height=11,width=7,unit="in",res=72
      )    
      draw(hm3, annotation_legend_list=lgd_list)
      dev.off()
  }
}