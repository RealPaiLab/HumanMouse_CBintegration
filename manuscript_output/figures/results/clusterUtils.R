library(patchwork)
library(ComplexHeatmap)
library(Seurat)
library(purrr)
library(tidyverse)
library(ggrepel)

# ' Run differential gene expression for all clusters.
run_all_de <- function(
  srat,
  clusters,
  grouping.var = "dataset_source",
  logfc.threshold = 0.2,
  min.pct=0.2
) {

# set Idents to UBC subclusters
Idents(srat) <- srat$seurat_clusters
  # run differential gene expression (each cluster vs. all other clusters)
  de_gene_list <- map(
    .x = clusters,
    .f = \(clust) {
      message(sprintf("***%s, n = %s cells***", clust, 
        sum(Idents(srat) == clust)))

      if (is.null(grouping.var)) {
        cat("No grouping variable provided, running DE without latent variables\n")
        # differential expression
        de_genes <- FindMarkers(
          object = srat,
          ident.1 = clust,
          assay = "SCT",
          logfc.threshold = logfc.threshold,
          min.pct = min.pct
        ) %>%
          # add columns for cluster and gene tested
          mutate(seurat_cluster = clust, .before = 1) %>%
          rownames_to_column(var = "gene")
      } else{
      # differential expression
      de_genes <- FindMarkers(
        object = srat,
        ident.1 = clust,
        grouping.var = grouping.var, 
        assay = "SCT",
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        test.use = "MAST",
        latent.vars = c(grouping.var)
      ) %>%
        # add columns for cluster and gene tested
        mutate(seurat_cluster = clust, .before = 1) %>%
        rownames_to_column(var = "gene")
      }

      # return differential expression results
      return(de_genes)
    }
  ) 

  de_genes <- do.call(rbind, de_gene_list)
  return(de_genes)
}

#' Plot complex heatmap of DE genes for clusters. Expects cluster information
#' to be stored in the "seurat_clusters" metadata column.
#'
#' @param srat Seurat object
#' @param de_genes Dataframe of differential expression results from run_all_de()
#' @param mdata Metadata dataframe with sample information
#' @param fcThresh (numeric) fold-change threshold for including genes in heatmap
#' @param QValThresh (numeric) Adjusted p-value threshold for including genes in heatmap
#' @param showTop (integer) number of top genes to label per cluster in heatmap
#' @param subsetGenes (character vector) optional vector of genes to subset DE results to before selecting top genes for heatmap. 
# This can be used to restrict the heatmap to e.g. protein-coding genes.
#' @param genesToInclude (character vector) genes to include in the heatmap regardless of their differential expression status
#' @return ComplexHeatmap object
plot_heatmap <- function(srat, de_genes,mdata, fcThresh=1.5, QValThresh=0.05, showTop=100, subsetGenes=NULL, genesToInclude=c()) {

  cat(sprintf("Subsetting DE genes to protein-coding genes\n"))
  if (!is.null(subsetGenes)) {
    de_genes <- de_genes %>%
      filter(gene %in% subsetGenes)
  }

  # rank top genes by log2FC
  feature_ranking <- list()
  for (clust in levels(srat$seurat_clusters)) {
    x <- de_genes %>%
      filter(seurat_cluster == clust, 
        p_val_adj < QValThresh, 
        avg_log2FC > log2(fcThresh)
      ) 
      feature_ranking[[clust]] <- x %>%
      mutate(rank = min_rank(desc(avg_log2FC))) %>%
      arrange(rank) %>%
      dplyr::select(gene, rank)
  }

  # get unique set of top genes across all clusters
unique_markers <- c()
  for (clust in names(feature_ranking)) {
    nr <- min(nrow(feature_ranking[[clust]]), showTop)
    feature_ranking[[clust]] <- feature_ranking[[clust]] %>%
      filter(rank <= nr)
    unique_markers <- c(unique_markers, feature_ranking[[clust]]$gene)
  }
unique_markers <- unique(unique_markers)

  # get scaled expression matrix for these genes
  xpr <- AverageExpression(
    object = srat,
    assays = "SCT",
    features = unique_markers,
    group.by = "seurat_clusters",
    slot = "data"
  )$SCT
  xpr <- t(scale(t(xpr)))

  # heatmap colours
  col_fun <- circlize::colorRamp2(
  breaks = c(-2,0,2),
  colors = c("#FF00FF","#000000", "#FFFF00")
  )

  # label top genes
  anno_genes <- map(
    .x = feature_ranking,
    .f = \(df) {slice_min(df, order_by = rank, n = 5)}
  ) %>%
    unlist() %>%
    unname() %>%
    c("OTX2", "CBFA2T2", "SOX4", "SOX11", "LMX1A","TOP2A","MKI67","DCX","NEUROD1")
  ha <- rowAnnotation(
    foo = anno_mark(at = match(anno_genes, rownames(xpr)), labels = anno_genes,
    labels_gp=gpar(fontsize = 24, col = "black"))
)
  cat("Top genes:\n")
  # print top genes per clusters
  for (clust in names(feature_ranking)) {
    top_genes <- feature_ranking[[clust]] %>%
      slice_head(n = 20)
    cat(sprintf("\tg%s: %s\n", clust, paste(top_genes$gene, collapse=", ")))
    cat("\n")
  }

  hm <- Heatmap(
    matrix = xpr,
    name = "Scaled Expression",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    col = col_fun,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 20),
    right_annotation = ha,
    column_title = "UBC Subclusters",
    row_title = "DE Genes" 
  )

  hm
}

# write a subroutine that takes as input a Seurat object and group column, 
# and produces DimPlots that highlight each group separately in dark red, 
# with the other groups in light gray
highlight_groups_dimplot <- function(srat, group_col) {
  groups <- sort(unique(srat[[group_col]][,1]))
  plotlist <- list()
  for (grp in groups) {
    #cat(sprintf("\tGroup: %s\n", grp))
    srat$highlight <- ifelse(srat[[group_col]][,1] == grp, grp, "Other")
    p <- DimPlot(srat, group.by = "highlight", cols = c("darkred", "lightgray"),order=TRUE) + 
      ggtitle(sprintf("CL %s", grp))

    p <- p + labs(x=NULL,y=NULL) + theme(legend.position="none")

      plotlist[[grp]] <- p
  }
  #p <- ggpubr::ggarrange(plotlist = plotlist, ncol = 4, nrow = ceiling(length(plotlist)/4))
   p <- cowplot::plot_grid(plotlist = plotlist, ncol = 4)
  return(p)
}


# characterizes clusters at multiple resolutions
# @param srat (Seurat) Seurat object with clustering already performed on SCT assay
# @param outDir (string) output directory for saved files
# @param clRes (vector) vector of clustering resolutions to analyze
# @return none (saves DimPlots and differential expression results to files)
clusterWorkup <- function(srat,outDir,clRes=seq(0.1,1,0.05),
  grouping.var=NULL){

  DefaultAssay(srat) <- "SCT"
  md <- srat[[]]
 testnames <- paste0("SCT_snn_res.",clRes)
  missing <- testnames[!(testnames %in% colnames(md))]
if (length(missing)>0){
   stop (sprintf("ERROR: missing clustering resolutions in metadata: %s.\nYou may need to generate them first.\n", paste(missing, collapse=",")));
  }

  for (rr in clRes){
     cat ("\n", rr );

     curDir <- sprintf("%s/Res_%s", outDir, rr);
     if (!dir.exists( curDir )) {
       dir.create( curDir, recursive = FALSE);
     } 
     
     res <- paste0("SCT_snn_res.",rr);res;
     srat$seurat_clusters <- srat[[res]];# so we can use as a common column 
     clF <- sprintf("%s/clusters.%s.png", curDir, res)
     p <- highlight_groups_dimplot(srat, "seurat_clusters")
     ggsave(p,file=clF,width=16, height=16, units="in", dpi=300)
  }
##outList <- list()
  cat("\nFinding markers for all clusters at all resolutions\n")
  #srat@misc$FindAllMarkerList <- list()
  for (rr in clRes){
    curDir <- sprintf("%s/Res_%s", outDir, rr);
    res <- paste0("SCT_snn_res.",rr);res;
    srat$seurat_clusters <- srat[[res]];# so we can use as a common column 
    de_genes <- run_all_de(
      srat = srat,
      clusters = levels(srat$seurat_clusters),
      grouping.var = grouping.var
    )
    #outList[[res]] <- de_genes
    outFile <- sprintf("%s/DEG_allClusters_%s.tsv", curDir, res)
    write.table(de_genes,
      file = outFile,
      sep = "\t",
      col = TRUE,
      row = TRUE,
      quote = FALSE
    )
  }
  #outList
}



  #' Plot a cell by gene heatmap for selected genes
  #' @param srat Seurat object
  #' @param selGenes Vector of gene names to include in heatmap
  #' @return ComplexHeatmap object
  plot_heatmap_cells <- function(srat, selGenes) {
    library(ComplexHeatmap)
    library(circlize)

    # get expression data for selected genes
    exprData <- GetAssayData(srat, assay="SCT", slot="data")
    
    selGenes2 <- selGenes[selGenes %in% rownames(exprData)]
    cat(sprintf("Found %s/%s selected genes in expression data\n", length(selGenes2), length(selGenes)))
    selGenes <- selGenes2
    exprData <- as.matrix(exprData[selGenes, ])
    browser()

    cat(sprintf("Subsetting expression data to %s selected genes\n", length(selGenes)))
    cat(sprintf("and %s cells\n", ncol(exprData)))
    
    if (nrow(exprData) == 0) {
      stop("No expression data found for selected genes")
    }

    # scale expression data by gene
    exprDataScaled <- t(scale(t(exprData)))

    # create annotation for clusters
    clusterAnno <- HeatmapAnnotation(
      df = data.frame(Cluster = srat$seurat_clusters),
      col = list(Cluster = structure(
        rainbow(length(levels(srat$seurat_clusters))),
        names = levels(srat$seurat_clusters)
      ))
    )

    # create heatmap
    hm <- Heatmap(
      exprDataScaled,
      name = "Scaled Expression",
      top_annotation = clusterAnno,
      show_row_names = TRUE,
      show_column_names = FALSE,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    )

    return(hm)
  }


#' Plot volcano plot of DE genes per cluster given the output of run_all_deg
#' @param deg (data.frame) of DE results from run_all_de()
#' @param showTopGenes (integer) how many top genes to highlight
#' @param logFCcutoff (numeric) the cutoff for visual depiction of upregulated or downregulated genes
#' @param titlePfx (string) prefix to add to the title of each plot
#' @param base_size (integer) base_size for theme() command
#' @return list of ggplot objects (volcano plot), one per cluster
plotVolcano_usingDEG <- function(deg, showTopGenes=50, logFCcutoff=0.25, titlePfx="", base_size=18) {
# plot a volcano plot for each cluster tested. label the top genes in each.
# color upregulated genes in red and downregulated genes in blue. use a significance threshold of 0.05 for adjusted p-value and a log2 fold change threshold of 0.25 for labeling genes.
pList <- list()
for (cl in unique(deg$seurat_cluster)) {
    deg_cl <- deg %>% filter(seurat_cluster == cl)
     deg_cl <- deg_cl %>%
      mutate(
        significant = case_when(
          p_val_adj < 0.05 & avg_log2FC > logFCcutoff ~ "upregulated",
          p_val_adj < 0.05 & avg_log2FC < -1*logFCcutoff ~ "downregulated",
          TRUE ~ "not_significant"
        )
      )
        p <- ggplot(deg_cl, aes(x=avg_log2FC, y=-log10(p_val_adj), color=significant)) +
      geom_point(size=0.5) +
      scale_color_manual(values = c("upregulated" = "red", 
        "downregulated" = "blue", "not_significant" = "grey")) +
      labs(x = "Average log2 fold change", y = "-log10 adjusted p-value") +
      ggtitle(sprintf("%s: Cluster %s", titlePfx, cl)) +
      theme_minimal(base_size = base_size) +
      theme(legend.position = "none")

    # label top 50 genes by adjusted p-value
    top_genes <- deg_cl %>% arrange(p_val_adj) %>% head(showTopGenes) %>% pull(gene)
    # use ggrepel to label the top genes on the volcano plot
    p <- p + geom_text_repel(data = deg_cl %>% filter(gene %in% top_genes), 
      aes(label = gene), size = 3, max.overlaps = Inf)
    #ggsave(sprintf("%s/volcano_human_cluster_%s.pdf", out_dir, cl), plot = p, width = 6, height = 5)
    pList[[paste0("cluster",cl)]] <- p
}
return(pList)
}
