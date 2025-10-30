

#' Runs DimPlot at all resolutions
#' @param sratIn (Seurat) Seurat object to plot
#' @param resSet (numeric) cluster resolutions to plot
#' @param origCluster (char) metadata column for colour-coding the 
#' original dimplot e.g., could be known cell types
#' @param outDir (char) output directory for files 
runDimPlot_multiRes <- function(sratIn, resSet, origCluster,outDir) {
    plotlist <- list()
    plt <- DimPlot(sratIn, reduction = "umap", label = TRUE, 
        repel = TRUE, group.by = origCluster)

    for (i in resSet) {
      out_folder <- paste0(outDir, "/", i, "_rl_res")
      print(out_folder)
      pattern <- paste0("snn_res.", i)
      if (!dir.exists(out_folder)) {
        dir.create(out_folder, recursive = FALSE)
      }
      dt <- format(Sys.Date(),"%y%m%d")

      #Plotting the UMAP of all clusters in resolution
      plt <- DimPlot(sratIn, reduction = "umap", 
        group.by = pattern, label = TRUE, label.size = 3)
      plt <- plt + ggtitle(paste("Resolution:", i))

      ggsave(filename = paste0("sratIn_", i, "_clustering_umap.png"), 
        plot = plt, path = out_folder, width = 7, 
        height = 7, units = "in", dpi = 300)
}

#' Plots DimPlot such that each cluster in a resolution is shown singly
#' highlighted in deep red, against a background of gray.
#' @param sratIn (Seurat) Seurat object to plot
#' @param resSet (numeric) cluster resolutions to plot
#' @param outDir (char) output directory for files
runDimPlot_singleClusterView_allRes <- function(sratIn, resSet, outDir) {
    # Plotting individual clusters
    for (curRes in resSet) {
        cat(sprintf("Resolution = %1.2f",curRes))
        curCol <- paste0("snn_res.", i)

        plist <- list()
        out_folder <- paste0(outDir, "/", i, "_rl_res")

        for (k in levels(sratIn@meta.data[[curCol]])) {
            sratIn$cur_cluster <- sratIn@meta.data[[curCol]] == k
            plist[[k]] <- DimPlot(sratIn  , reduction="umap", 
                cols=c("grey90","#aa0000"),
                group.by="cur_cluster",
                order=TRUE,
                label=FALSE) + ggtitle(k) + 
                xlab("UMAP 1") + ylab("UMAP 2") + 
                theme(plot.title=element_text(size=9)) + NoLegend()
        }
        p <- wrap_plots(plist, ncol = 3)
        ggsave(filename = paste0(curRes, "_individual_clusters_umap.png"), 
            plot = p, path = out_folder, width = 9, 
            height = ceiling(length(levels(sratIn@meta.data[[pattern]]))/3)*3, 
            units = "in", dpi = 300
        )
    }
}

#' Plots DimPlot such that each cluster in a resolution is shown singly
#' highlighted in deep red, against a background of gray.
#' @param sratIn (Seurat) Seurat object to plot
#' @param resSet (numeric) cluster resolutions to plot
#' @param outDir (char) output directory for files
#' @param geneSet (char) vector of genes to plot
runDotPlot <- function(sratIn, resSet, outDir,geneSet){

}


runDEG <- function(){

}


#' plot heatmap of markers for each cell cluster.
#' @param srat (Seurat) Seurat object to plot heatmap
#' @param DEGresults (data.frame) table with FindMarkers() result combined
#' for all one-versus-all cluster comparisons
#' (e.g., cluster 0 vs all, cluster 1 vs all etc.,). Table must have a 
#' "avgLogFC" column with the average log fold change of the gene. This value
#' will be used for the cutoff
#' @param cutoffs (numeric) avgLog2FC cutoffs to call "markers"
#' @param cluster_order (char) level order for clusters. If NULL, the 
#' default order in the metadata column will be used
#' @return (list) Heatmap objects
#'  1. 

plotHeatmap <- function(srat, DEGresults,cutoffs=c(1,1.5,2),
    cluster_order=NULL) {    
    library(ComplexHeatmap)

    if (is.null(cluster_order)) {
        cat("cluster_order not provided. Sorting alphanumerically...\n")
        cluster_order <- sort(unique(DEGresults$cluster))
    }

    if (!("avgLogFC" %in% colnames(DEGresults))) {
        stop("DEGresults must have an 'avgLogFC' column for cutoffs. Please add one.")
    }

    for (cutoff in cutoffs) {
        message("---------------------------------")
        message(sprintf("cutoff = %1.2f", cutoff))

        topGenes <- list()
        topData <- list()
        message("collecting top  genes per UBC cluster")
        for  (cl in cluster_order) {         
            dat <- subset(orig_dat,cluster %in% cl)
            browser()
            dat <- dat[order(dat$avgLogFC,decreasing=TRUE),]
           # message(cl)
            #print(head(dat[,c("avgLogFC","feature")]))
             dat <- subset(dat,
                avgLogFC > cutoff)
            nr <- min(100,nrow(dat))
            if (nr > 0) {
                topGenes[[cl]] <- dat$feature[1:nr]
                topData[[cl]] <- cbind(rep(cutoff,nr),rep(cl,nr),dat[1:nr,])
            }

        }
        topData <- do.call("rbind",topData)
        colnames(topData)[1:2] <- c("avgLogFC_cutoff","topGene_cluster")
        write.table(topData,
            file=sprintf("%s/topGenes.cutoff%1.2f.%s.txt",
                outDir,cutoff,fileSfx),
            sep="\t",col=TRUE,row=FALSE,quote=FALSE
        )

        message("\tplotting heatmap at cellular resolution")
        feat <- c(unlist(topGenes),
                "SOX4","SOX11","FOXP2","EYS","CRX","NRL")

        hm <- DoHeatmap(
            srat_ubc,
            features = feat,
            group.by = "broad_w_ubc_subtypes",
            size=3
        ) + theme(axis.text.y=element_blank())

        hm <- hm + ggtitle(sprintf("UBC markers (Log2FC > %1.2f)",cutoff))
        ggsave(hm, file=sprintf("%s/heatmap_UBCmarkers_cutoff%1.2f.%s.png",
            outDir,cutoff, dt)
        )

        message("plot heatmap at cluster resolution")
        xpr <- as.matrix(AverageExpression(srat_ubc,
            features=feat,
            assay="RNA",
            group.by="broad_w_ubc_subtypes"
        )[[1]])
        colnames(xpr) <- sub("-","_",colnames(xpr))
        xpr <- xpr[,cluster_order]

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
                column_title=sprintf("UBC markers (log2FC > %1.2f), unscaled",
                    cutoff), 
            top_annotation = topannot,
            right_annotation = ha
        )

        png(sprintf("%s/heatmap_AveUBCmarkers_cutoff%1.2f.unscaled.%s.png",
            outDir, cutoff,dt),
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
                column_title=sprintf("UBC markers (log2FC > %1.2f), scaled",
                    cutoff),
                show_row_names = FALSE,
                top_annotation = topannot,
                right_annotation = ha
        )

        png(sprintf("%s/heatmap_AveUBCmarkers_cutoff%1.2f.scaled.%s.png",
            outDir, cutoff,dt),
            height=11,width=7,unit="in",res=72
        )    
        draw(hm3, annotation_legend_list=lgd_list)
        dev.off()
    }
}