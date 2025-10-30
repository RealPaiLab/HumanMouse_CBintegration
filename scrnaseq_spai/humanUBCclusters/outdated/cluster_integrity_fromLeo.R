
library(tidyverse)
library(Seurat)
library(qs)
library(RColorBrewer)
library(patchwork)
library(ggrepel)
library(doParallel)
library(foreach)
library(ComplexHeatmap)
library(ggpubr)

#' @param dataset_srat (Seurat) dataset to process
#' @param outDir (char) path to directory to save all output logs & graphs
#' @param runSCT (logical) run SCTransform?
#' @param runDEG (logical) run FindMarkers() on all clusters?
#' @param plotClusters (logical) Run DimPlots?
#' @param resSet (numeric) vectors of resolutions for FindClusters()
#' @param dotplot_genes2plot (char) genes for DotPlot
#' @param featureplot_genes2plot (char) genes for FeaturePlot
#' @param origClusterColumn (char) metadata column to colour code the DimPlot
#' @param heatmap_topAnnot_columns (char) metadata columns to show at the top
#' of the heatmap with cluster markers, showing breakdown. 
#' e.g., age, sex distribution in the cluster
#' at a given resolution (e.g., putative cell assignments)
#' @param heatmap_log2FC_cutoffs (numeric) log2FC cutoffs for DEG analysis to 
#' plot the heatmap
#' @param heatmap_addGenes (char) additional genes to show on DEG heatmap
#' @param gencodeList (data.frame) gene definitions. If provided, will filter
#' markers for just these.
profileClusters <- function(dataset_srat, outDir, 
    runSCT=FALSE, runDEG=TRUE, 
    resSet=seq(0.2,0.8,0.2),dotplot_genes2plot, origClusterColumn,
    heatmap_topAnnot_columns, heatmap_log2FC_cutoffs=c(1,1.5,2),
    heatmap_addGenes=c("SOX4","SOX11","MKI67","CRX","NRL"), 
    featureplot_genes2plot,
    numCores=4L,gencodeList=NULL) {

    set.seed(1)
    source("cluster_multiple_res.R")
    out_directory <- outDir

    cat(sprintf("About to profile clusters = { %s }\n",
        paste(resSet, collapse=","))
    )

    if (runSCT){
    cat("\n\n\n")
    cat("-------------------------------------\n")
    cat(sprintf("Running SCTransform\n"))
    cat("-------------------------------------\n")
    dataset_srat <- SCTransform(dataset_srat, 
        vars.to.regress = c("CC.Difference", "sex"), 
        variable.features.n = 5000)
    
    } else {
        cat("** Skipping SCTransform.\n")
    }

    cat("\n\n\n")
    cat("-------------------------------------\n")
    cat(sprintf("Running PCA\n"))
    cat("-------------------------------------\n")
    dataset_srat <- RunPCA(dataset_srat)

    cat("-------------------------------------\n")
    cat(sprintf("Running ScaleData() ...\n"))
    cat("-------------------------------------\n")
    dataset_srat <- ScaleData(dataset_srat,
        assay="SCT",features=rownames(dataset_srat)
    )

    cat("\n\n\n")
    cat("-------------------------------------\n")
    cat(sprintf("Running ElbowPlot() ...\n"))
    cat("-------------------------------------\n")
    plt <- ElbowPlot(dataset_srat, ndims = 50, reduction = "pca") +            theme_classic()
    
    ggsave(filename = "dataset_srat_elbow_plot.png", 
      plot = plt, path = out_directory, width = 7, 
      height = 7, units = "in", dpi = 600)

    cat("\n\n\n")
    cat("-------------------------------------\n")
    cat(sprintf("Running UMAP with 25 dims ...\n"))
    cat("-------------------------------------\n")

    dataset_srat <- RunUMAP(dataset_srat, dims = 1:25) %>% 
      FindNeighbors(dims = 1:25)

    # Removing previously generated clusters
    pattern <- paste0(".*_snn_res.")
    matching_column <- grep(pattern, 
        colnames(dataset_srat@meta.data), value = TRUE)
    for(column in matching_column){
        dataset_srat[[column]] <- NULL
    }

    cat("\n\n\n")
    cat("-------------------------------------\n")
    cat(sprintf("Clustering at multiple resolutions ...\n"))
    cat("-------------------------------------\n")

    # The last two parameters for this function is just denotating the  integration type and the out directory which we don't need to use here
    resolutions <- resSet
    dataset_srat <- cluster_multiple_res(
      dataset_srat,
      resolutions,
      "",
      ""
    )

    cat("\n\n\n")
    cat("-------------------------------------\n")
    cat(sprintf("Plotting Dim Plot of original cluster assignments ...\n"))
    cat("-------------------------------------\n")
    plt <- DimPlot(dataset_srat, reduction = "umap", label = TRUE, 
        repel = TRUE, group.by = origClusterColumn)

    ggsave(filename = sprintf("%s_umap.png", origClusterColumn),
      plot = plt, path = out_directory, width = 7, 
      height = 7, units = "in", dpi = 600)

    cat("\n\n\n")
    cat("-------------------------------------\n")
    cat(sprintf("Plotting FeaturePlot of original cluster assignments ...\n"))
    cat("-------------------------------------\n")
    p <- FeaturePlot(dataset_srat, featureplot_genes2plot,order=TRUE, ncol=3)
    ggsave(p, file=sprintf("%s/FeaturePlot_%s.png",out_directory,dt),
          width = 9, 
          height = ceiling(length(featureplot_genes2plot)/3)*3, units = "in")

    cat("\n\n\n")
    cat("-------------------------------------\n")
    cat(sprintf("Plotting DimPlot and DotPlot at different resolutions ...\n"))
    cat("-------------------------------------\n")

    #cl <- makeCluster(numCores)
    #registerDoParallel(numCores)

    for (i in resolutions) { #%do% {
      cat(sprintf("\nResolution = %1.2f\n",i))
      out_folder <- paste0(out_directory, "/", i, "_res")
      print(out_folder)
      pattern <- paste0("snn_res.", i)
      if (!dir.exists(out_folder)) {
        dir.create(out_folder, recursive = FALSE)
      }
      dt <- format(Sys.Date(),"%y%m%d")

      cat("\t* Plotting the UMAP of all clusters in resolution\n")
      plt <- DimPlot(dataset_srat, reduction = "umap", group.by = pattern, label = TRUE, label.size = 3) + ggtitle(paste("Resolution:", i))

      ggsave(filename = paste0("res", i, "_clustering_umap.png"), 
      plot = plt, path = out_folder, width = 7, 
      height = 7, units = "in", dpi = 600)

      cat("\t* Plotting clusters singly\n")
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
        ggsave(filename = paste0("res", i, "_individual_clusters_umap.png"), 
          plot = p, path = out_folder, width = 9, 
          height = ceiling(length(levels(dataset_srat@meta.data[[pattern]]))/3)*3, units = "in", dpi = 600)
    
      cat("\t* Plotting DotPlot\n")
      # Plotting dot plot
      #genes <- c("MKI67","WLS", "EOMES", "LMX1A", "OTX2", "RBFOX3", "ATOH1","PAX6", "RELN")
      Idents(dataset_srat) <- pattern

      # Note: Plotting the SCT assay
      plt <- DotPlot(object = dataset_srat, assay = "SCT", 
        features = dotplot_genes2plot) + 
        theme_classic() 
      ggsave(filename = paste0("res", i, "_dot_plot.png"), 
        plot = plt, path = out_folder, width = 10, 
        height = 5, units = "in", dpi = 600)

     cat("\t* Running differential expression\n")
      # Differential expression

      # Get counts slot from Seurat object
      counts <- GetAssayData(object = dataset_srat, slot = "counts")


     DEGresultsFile <- paste0(out_folder, "/", i, "_markers.csv")
     if (runDEG){
        # If a count for a gene in a cell is greater than 0, set as TRUE (= 1)
        nonzero <- counts > 0
        # If 1% or more cells are TRUE, keep the gene. Each TRUE value = 1. Taking the sum of all the cells for that gene
        keep_genes <- Matrix::rowSums(nonzero >= (0.01*length(Cells(dataset_srat))))

        # Get the genes names that we are keeping
        keep_genes <- rownames(dataset_srat)[keep_genes]

        ## # Note that the number of genes in the Aldinger and Sepp datasets are the same as they were merged together before
        cat(sprintf("\t\t Filtering for genes expressed in >1%% cells: %i total, %i in 1%%\n", 
            length(rownames(dataset_srat)),length(keep_genes)))
        dataset_srat <- subset(dataset_srat, features = keep_genes)

        cat("\t\tPrepSCTFindMarkers\n")
        dataset_srat <- PrepSCTFindMarkers(object = dataset_srat)
        cat("\t\tReleveling\n")
        # Releveling based on alphabetical then numeric order
        levels_index <- str_order(
          unique(dataset_srat@meta.data[[pattern]]),
          numeric = TRUE
        )
        dataset_srat@meta.data[[pattern]] <- fct_relevel(dataset_srat@meta.data[[pattern]], levels(dataset_srat@meta.data[[pattern]][levels_index]))

        # Setting the Idents to be the cluster name because FindMarkers() uses the Ident
        Idents(dataset_srat) <- pattern

        cat("\t\tRunning DEG analysis\n")
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
              assay = "SCT",
              recorrect_umi = FALSE
            )
          # Append markers to list
          marker_list[[cluster]] <- markers
        }
        cat("\t\tmerging in DF\n")
        # Creating a dataframe from the list of markers
        df_cells <- marker_list
        for (name in names(df_cells)) {
          # Adding cluster name as a column
          df_cells[[name]]$cluster <- name
          # Adding gene name as a column
          df_cells[[name]]$feature <- rownames(df_cells[[name]])
        }
        DEGresults <- do.call(rbind, df_cells)
        write.csv(DEGresults, 
          file = DEGresultsFile, row.names = FALSE)
        cat("\t\tDEG analysis done. Results saved.\n")
     } else {
        cat("\t\t** Skipping DEG FindMarker() calling. Looking for results file...")
        DEGresults <- read.delim(DEGresultsFile,sep=",",h=TRUE)
     }
      levels_index <- str_order(
          unique(dataset_srat@meta.data[[pattern]]),
          numeric = TRUE
      )
      
      DEGresults$cluster <- as.integer(DEGresults$cluster)

       DEGresults$ptrans <- -log10(DEGresults$p_val)
       maxp <- max(DEGresults$ptrans[which(!is.infinite(DEGresults$ptrans))])
       DEGresults$ptrans[which(is.infinite(DEGresults$ptrans))] <- maxp

      cat(sprintf("\t\t Plotting volcano plot (all genes)"))
      plist <- list()
      plist2 <- list()
      for (curCluster in unique(DEGresults$cluster)){
        curDEG <- subset(DEGresults, cluster %in% curCluster)
       
        p <- ggplot(curDEG,aes(x=avg_log2FC,y=ptrans))
        p <- p + geom_point(aes(colour=avg_log2FC),size=0.5) +
            scale_colour_gradient2()
        p <- p + xlab("-log10(p-value)") + ylab("log2FC") 

        idx <- is.infinite(curDEG$ptrans)
        if (any(idx)) curDEG[idx] <- 300
        p <- p + xlim(min(curDEG$avg_log2FC)*1.15, max(curDEG$avg_log2FC)*1.15)
        p <- p + ylim(0, max(curDEG$ptrans)*1.15)
        p <- p + ggtitle(sprintf("%s: Cluster %i",pattern, curCluster))

        curDEG$absLogFC <- abs(curDEG$avg_log2FC)
        curDEG <- curDEG[order(curDEG$absLogFC,decreasing=TRUE),]
        curDEG2 <- subset(curDEG, p_val_adj < 0.05 & absLogFC > 2)
        if (nrow(curDEG2)>1) {
        p <- p + geom_label_repel(
           data=curDEG2, 
          aes(x=avg_log2FC, y=ptrans, label=feature),
          max.overlaps=100, 
          size=1) #+ xlim(c(0,1.5)) + ylim(c(0,2))
        }

        plist[[as.character(curCluster)]] <- p

        # GENCODE genes
        curDEG <- subset(curDEG, feature %in% gencodeList$gene_name)
         p <- ggplot(curDEG,aes(x=avg_log2FC,y=ptrans))
        p <- p + geom_point(aes(colour=avg_log2FC),size=0.5) +
            scale_colour_gradient2()
         p <- p + xlim(min(curDEG$avg_log2FC)*1.15, max(curDEG$avg_log2FC)*1.15)
        ymax <- min()
       p <- p + ylim(0, max(curDEG$ptrans)*1.15)
        p <- p + xlab("avg log2FC") + ylab("-log10(p)") 
        p <- p + ggtitle(sprintf("%s: Cluster %i",pattern, curCluster))

        curDEG$absLogFC <- abs(curDEG$avg_log2FC)
        curDEG <- curDEG[order(curDEG$ptrans,decreasing=TRUE),]
        curDEG2 <- subset(curDEG, p_val_adj < 0.05 & absLogFC > 1)
        if (nrow(curDEG2)>1) {
          nr <- min(35, nrow(curDEG2))
          print(nr)
          curDEG2 <- curDEG2[1:nr,]
        p <- p + geom_text_repel(
           data=curDEG2, 
          aes(x=avg_log2FC, y=ptrans, label=feature),
          max.overlaps=100, 
          size=3) #+ xlim(c(0,1.5)) + ylim(c(0,2))
        
        plist2[[as.character(curCluster)]] <- p
        }   
      }
      x <- ggarrange(plotlist=plist)
      ggsave(x, file=sprintf("%s/volcano.%s.png",
              out_folder, dt),
              width = 12, 
            height = ceiling(length(levels(dataset_srat@meta.data[[pattern]]))/3)*4, 
            units = "in", dpi = 600)

      x <- ggarrange(plotlist=plist2)
      ggsave(x, file=sprintf("%s/volcano.Gencodeonly.%s.png",
              out_folder, dt),
              width = 12, 
            height = ceiling(length(levels(dataset_srat@meta.data[[pattern]]))/3)*4, 
            units = "in", dpi = 600)
  


      if (!is.null(gencodeList)) {
        cat("\t\t**** Filtering for GENCODE definitions ****\n")
        DEGresults <- subset(DEGresults, feature %in% gencodeList$gene_name)
        nr <- nrow(DEGresults)
        cat(sprintf("\t\t*** DEG file = %i entries; after filtering = %i genes\n *** ",
            nr, nrow(DEGresults)))
      }
    
      cat("\t\tPlotting heatmap -- TOFIX: make annotation general.\n")
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

      blah <- dataset_srat@meta.data[,c(pattern,"age")]
      blah$numage <- as.numeric(sub(" PCW","",blah$age))
      x <- aggregate(blah$numage, by=list(blah[,pattern]),FUN=mean,na.rm=TRUE)
      levels_index <- order(x[,2]) # order by increasing average age

      # reorder all
      # Releveling based on alphabetical then numeric order
      ages <- ages[levels_index,]
      sex <- sex[levels_index,]
      clusterOrder <- x[levels_index,1]

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
    
      cat(sprintf("\t\tPlotting heatmap for different cutoffs\n"))
      cutoffs <- heatmap_log2FC_cutoffs
      for (cutoff in cutoffs) {
          cat("---------------------------------")
          cat(sprintf("\t\t\tcutoff = %1.2f", cutoff))

          topGenes <- list()
          topData <- list()
          cat("\t\t\tcollecting top genes per cluster")
          
          
          for  (cl in clusterOrder) {
              dat <- subset(DEGresults,cluster %in% cl)
              dat <- subset(dat, 
                  p_val_adj < 0.05 
              )
              # Taking p_val, avg_log2FC, and feature columns
              dat <- dat[order(dat$avg_log2FC,decreasing=TRUE),]
            # cat(cl)
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
              file=sprintf("%s/TopGenes.cutoff%1.2f.%s.txt",
                  out_folder,cutoff,dt),
              sep="\t",col=TRUE,row=FALSE,quote=FALSE
          )

          cat("\tplotting heatmap at cellular resolution")
          feat <- c(unlist(topGenes), heatmap_addGenes)

          hm <- suppressWarnings(DoHeatmap(
              dataset_srat,
              features = feat,
              group.by = pattern,
              size=3
          )) + theme(axis.text.y=element_blank())

          hm <- hm + ggtitle(sprintf(
                "Top cluster genes (Log2FC > %1.2f)",cutoff))
          ggsave(hm, file=sprintf("%s/heatmap_TopGenes_cutoff%1.2f.%s.png",
              out_folder,cutoff, dt)
          )

          cat("plot heatmap at cluster resolution")
          xpr <- as.matrix(AverageExpression(
                dataset_srat,
                features=feat,
                assay="SCT",
                group.by=pattern
          )[[1]])
          colnames(xpr) <- sub("-","_",colnames(xpr))
          xpr <- xpr[,levels_index]

          quants <- quantile(xpr,c(0.1,0.95))
          col_fun <- circlize::colorRamp2(quants,c("#000000","#FFFF00"))

          annoGenes <- heatmap_addGenes
          cat("Genes to annotate")
          for (nm in 1:length(topGenes)) {
              cur <- topGenes[[nm]]; ln <- length(cur)
              cur <- cur[1:min(length(cur),5)]
              annoGenes <- c(annoGenes, cur)
              cat(sprintf("%s: %i in set; selected = { %s }", 
                  nm, ln, paste(cur,collapse=",")))
          }

          annoGenes <- c("MKI67","WLS","EYS","IQCJ-SCHIP1","EOMES","LRB1B","SYNE1","KIRREL3")

          # add top genes per cluster on the right
          ha <- rowAnnotation(foo = 
            anno_mark(at = match(annoGenes,rownames(xpr)),
              labels = annoGenes))
          hm2 <- Heatmap(xpr, name="AverageExpression",
                  cluster_rows=FALSE,
                  cluster_columns=FALSE,
                  col = col_fun,
                  show_row_names = TRUE,
                  column_title=sprintf("Top genes (log2FC > %1.2f), unscaled",
                      cutoff), 
              top_annotation = topannot,#
              right_annotation = ha
          )

          png(sprintf("%s/heatmap_AvgXprTopGenes_cutoff%1.2f.unscaled.%s.png",
              out_folder, cutoff,dt),
              height=11,width=7,unit="in",res=72
          )
          draw(hm2, annotation_legend_list = lgd_list)
          dev.off()

          cat("now plot scaled")
          unscaled <- xpr
          xpr <- t(scale(t(xpr)))
          col_fun <- circlize::colorRamp2(
              c(-2,0,2),c("#FF00FF","#000000", "#FFFF00")
            #c(0, max(xpr)), c("#000000","#FFFF00")
          )

          # add top genes per cluster on the right
          ha <- rowAnnotation(foo = 
            anno_mark(at = match(annoGenes,rownames(xpr)),
              labels = annoGenes))
          hm3 <- Heatmap(xpr, name="AverageExpression",
                  cluster_rows=FALSE,
                  cluster_columns=FALSE,
                  col = col_fun,
                  column_title=
                    sprintf("Top genes (log2FC > %1.2f), scaled",
                      cutoff),
                  show_row_names = TRUE,
                  row_names_gp = gpar(fontsize=6),
                  top_annotation = topannot, #@ ,
                  right_annotation = ha
          )

          png(sprintf("%s/heatmap_AvgExpr_cutoff%1.2f.scaled.%s.png",
              out_folder, cutoff,dt),
              height=11*2,width=7*2,unit="in",res=300
          )    
          draw(hm3, annotation_legend_list=lgd_list)
          dev.off()
          cat("finished drawing all heatmaps.")

          cat("\n\n")
      }
    }
}
