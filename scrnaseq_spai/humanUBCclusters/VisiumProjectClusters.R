# project UBC clusters onto Visium data
rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggcorrplot)

# devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
#library(spacexr)

UBC_Seurat <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"
visiumFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Aldinger_Visium_2021Data/obj.list.rds"

DEGdir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/diffExpr"
useDEfile <- sprintf("%s/250424/de_genes.tsv", DEGdir)
useFRfile <- sprintf("%s/250424/cluster_marker_ranking.rds", DEGdir)
clusterNames <- "UBC_clusterNames.txt"
heatmap_logFC_thresh <- 1.5 # logFC threshold for heatmap
heatmap_pval_thresh <- 0.05 # adjusted p-value threshold for heatmap

de_genes <- read.delim(file = useDEfile,sep="\t",h=T,as.is=T)
feature_ranking <- readRDS(file = useFRfile)
ubc <- read.delim(file = clusterNames, sep="\t", h=T, as.is=T)

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/HumanUBC/Aldinger_Visium"
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = T)
}

logFile <- sprintf("%s/log.txt", outDir)
sink(logFile, split=TRUE)

tryCatch({
rotate_image <- function(p, rot_angle) {
    gt <- ggplot_gtable(ggplot_build(p))
    panel_idx <- which(gt$layout$name == "panel")
    rot_vp <- viewport(angle = rot_angle)
    gt[["grobs"]][[panel_idx]] <- editGrob(gt[["grobs"]][[panel_idx]], vp = rot_vp)
    p_rot <- ggdraw() + draw_grob(gt)

          return(p_rot)
      }
          #cat("Reading UBC clusters ...")
    #cat(sprintf("Input file:\n%s\n", UBC_Seurat))
    #t0 <- Sys.time()
    #srat <- readRDS(UBC_Seurat);
    #cat("done in ", Sys.time()-t0, "\n")

    ## data prep ##
    cat("Reading Visium data\n")
    
    cat(sprintf("Visium file:\n%s\n", visiumFile))

    visium <- readRDS(visiumFile)
        lbls <- c("sample1_slice1", "sample1_slice2", "sample2_slice1", "sample2_slice2")
    useSlices <- 1:4
    rotateAngles <- c(90, 90, -90, -90)

    cat("Making UBC subfolders\n")
    for (ubcGroup in ubc$cluster) {
      dir.create(sprintf("%s/%s", outDir, ubc$name[which(ubc$cluster == ubcGroup)]), showWarnings = FALSE)
    }


    cat("Plotting marker genes\n")
    # check slices
    geneList <- list(
      WLS="RLVZ",    
      EOMES="UBC",
      RORA="Purkinje cell",
      PAX6="GCP",
      GDF10="Bergmann Cell",      
      PAX2="GABAergic",  
      OLIG2="OPC",
      AQP4="Astrocyte",
      RBFOX3="Granule cell",
      MKI67="progenitor"
    )

cat("Plotting celltype markers")
   for (k in useSlices){
     #   cat(sprintf("Processing slice %s\n", lbls[k]))
        plotlist <- list()
        geneList2 <- intersect(names(geneList), rownames(visium[[k]]))
        for (gene in geneList2){
          short_lbl <- sub("sample", "S", lbls[k])
          short_lbl <- sub("_slice", "sl", short_lbl)
         ## cat(sprintf("Plotting %s\n", gene))
          p <- SpatialFeaturePlot(
            visium[[k]], features = gene)
          p <- p + ggtitle(sprintf("%s: %s (%s)", 
            short_lbl, gene, geneList[[gene]])) +
            theme(plot.margin = unit(c(0,0,0,0), "cm"),
              legend.position = "none") 
  
          plotlist[[gene]] <-rotate_image(p,rotateAngles[k])
          #ggsave(sprintf("%s/%s_%s.png", outDir, lbls[k],gene), 
          #p, dpi = 300)
          ##cat("saved\n")
       # visium[[k]]@meta.data$orig.ident <- lbls[k]
      }
      outFile <- sprintf("%s/Visium_%s_cellTypeGenes.png", 
      outDir, lbls[k])
      cat(sprintf("Saving to %s\n", outFile))
      x <- patchwork::wrap_plots(plotlist,ncol=5)
      suppressMessages(ggsave(x, file=outFile, width=18, height=12))
    }

cat("Now computing the correlation between EOMES and UBC markers\n")

topGenes <- 25
pvalCutoff <- 0.05/(topGenes * length(ubc$cluster)) # Bonferroni correction
corrCutoff <- 0.01

EOMESsigGenes <- list() # genes that significantly correlate with EOMES
for (ubcGroup in ubc$cluster) {
      cat("\nLooking at UBC cluster ", ubcGroup, "\n")
      ubcName <- ubc$name[which(ubc$cluster == ubcGroup)]

      cur <- feature_ranking[[ubcGroup]]
      geneList <- c(head(cur$gene,topGenes),"EOMES")
      if (ubcGroup == "UBC_5") {
        geneList <- c(geneList, "RFX3")
      } else if (ubcGroup == "UBC_8"){
        geneList <- c(geneList, "SMAD9")
      } else if (ubcGroup == "UBC_4"){
        geneList <- c(geneList, "TACR1")
      }

      cat(sprintf("%s: Tested genes: %s\n", 
        ubcName, paste(geneList, collapse=", ")))

      EOMESsigGenes[[ubcGroup]] <- list()
      for (k in 1:4) EOMESsigGenes[[ubcGroup]][[k]] <- c()
 
    #  "SOX2","WLS") #("MKI67","SOX2","ITPR1","EYS","PEX5L",) #("BRCA1")
      for (k in useSlices) {
        cat(sprintf("\n\nProcessing slice %s\n", lbls[k]))
        plotlist <- list()
        xpr <- as.matrix(visium[[k]][["SCT"]]@counts)
        geneList2 <- intersect(geneList, rownames(xpr))
        xpr <- t(xpr[geneList2,])

        corMat <- cor(xpr, method="pearson")
        corMat[is.na(corMat)] <- 0
      
        outFile <- sprintf("%s/%s/Visium_%s_%s_corrplot.pdf", 
          outDir, ubcName, ubcGroup, lbls[k])
        cat(sprintf("Saving to %s\n", outFile))
        testRes <- cor_pmat(xpr, method="pearson")
        testRes[is.na(testRes)] <- 1

         # which genes significantly correlate with EOMES?
        tmp <- 
          rownames(corMat)[which(testRes[,"EOMES"] < pvalCutoff &
          corMat[,"EOMES"] > corrCutoff)]
        tmp <- 
          setdiff(tmp, "EOMES")
        tmp <- corMat[tmp,"EOMES"]
        EOMESsigGenes[[ubcGroup]][[k]] <- sort(tmp, decreasing = TRUE)
      

        p <- ggcorrplot(corMat, p.mat = testRes, hc.order = TRUE,
          lab = TRUE, sig.level = pvalCutoff, # Bonferroni correction
          type="lower", insig = "blank",
          colors = c("red", "white", "blue"))
      
        p <- p + ggtitle(sprintf("Correlation plot: %s (%s)", 
            ubcName, lbls[k]))
        ggsave(filename = outFile, plot = p, width = 10,height =10 )
        }
  }
  cat("Genes that significantly correlate with EOMES:\n")
  for (ubcGroup in names(EOMESsigGenes)) {
    cat(sprintf("UBC cluster %s (%s): %s\n", 
      ubcGroup, 
      ubc$name[which(ubc$cluster == ubcGroup)],
      paste(names(EOMESsigGenes[[ubcGroup]]), collapse=", ")))
  }

consistentGenes <- list()
cat("Tallying which genes consistently correlate with EOMES\n")
for (k in names(EOMESsigGenes)){
  x <- unlist(EOMESsigGenes[[k]])
  nm <- table(names(x))
  nm <- sort(nm, decreasing = TRUE)
  nm <- nm[nm > 1]
  cat(sprintf("%s: Genes that show up in more than one slice\n",k))
  print(nm)
  consistentGenes[[k]] <- nm
}


cat("Plotting UBC cluster markers\n")
for (ubcGroup in ubc$cluster) {
      cat("Looking at UBC cluster ", ubcGroup, "\n")
      ubcName <- ubc$name[which(ubc$cluster == ubcGroup)]
      
    #  "SOX2","WLS") #("MKI67","SOX2","ITPR1","EYS","PEX5L",) #("BRCA1")
      for (k in useSlices){
        cat(sprintf("Processing slice %s\n", lbls[k]))

      geneSet <- consistentGenes[[ubcGroup]] #EOMESsigGenes[[ubcGroup]][[k]]
      geneList <- names(geneSet)
      if (length(geneList) < 1) {
        cat("No genes significantly correlated with EOMES. Skipping\n")
        next
      } else{
        geneList <- c("EOMES", geneList)
      }

        plotlist <- list()
        geneList2 <- intersect(geneList, rownames(visium[[k]]))
        for (gene in geneList2){
          short_lbl <- sub("sample", "S", lbls[k])
          short_lbl <- sub("_slice", "sl", short_lbl)
        #  cat(sprintf("Plotting %s\n", gene))
          p <- SpatialFeaturePlot(
            visium[[k]], features = gene)
          p <- p + ggtitle(sprintf("%s\n%s: %s (cor=%1.2f)", #(%s)", 
            ubc$name[which(ubc$cluster == ubcGroup)], short_lbl, gene,
            EOMESsigGenes[[ubcGroup]][[k]][which(names(geneSet) == gene)] #, geneList[[gene]]
            )) +
            theme(plot.margin = unit(c(0,0,0,0), "cm"),
              legend.position = "none")
          plotlist[[gene]] <- rotate_image(p,rotateAngles[k])    
       #   cat("saved\n")       
      }      
      outFile <- sprintf("%s/%s/Visium_%s_%s_ConsistentCorrGenes.png", 
      outDir, ubcName, ubcGroup, lbls[k])
      #cat(sprintf("Saving to %s\n", outFile))
      x <- patchwork::wrap_plots(plotlist,ncol=4)
      ggsave(x, file=outFile, width=12, height=12)
    }
}



####' deconvolute hindbrain visium against scRNA ref using RCTD
####' @param visium (SeuratObject) visium spatial object
#######' @param ref (SeuratObject) scRNA reference
#######' @param ref_cell_type (character) Column to be used as reference cell type
#######' @param max_cores (numeric) Max number of cores for RCTD
#######' @return (SeuratObject) The original visium object with deconvoluted cell ###type weigths
######deconv_hb_visium_RCTD <- function(visium, ###
######                                  ref, ref_cell_type, ###
######                                  max_cores = 8) {###
######  # extract information to pass to the RCTD Reference function###
#########
######  counts <- ref@assays$RNA$counts###
######  cluster <- as.factor(ref@meta.data[[ref_cell_type]])###
######  names(cluster) <- rownames(ref@meta.data)###
######  nUMI <- ref@meta.data$nCount_RNA###
######  names(nUMI) <- rownames(ref@meta.data)###
######  reference <- Reference(counts, cluster, nUMI)###
######  ###
######  # set up query with the RCTD function SpatialRNA###
######  counts <- visium[["Spatial"]]$counts###
######  coords <- GetTissueCoordinates(visium)###
######  colnames(coords) <- c("x", "y")###
######  coords[is.na(colnames(coords))] <- NULL###
######  query <- SpatialRNA(coords, counts, colSums(counts))###
######  ###
######  # RCTD###
######  RCTD <- create.RCTD(query, reference, max_cores = max_cores)###
######  RCTD <- run.RCTD(RCTD, doublet_mode = "full")###
######  ###
######  # assign cell type weights###
######  weights <- RCTD@results$weights###
######  norm_weights <- normalize_weights(weights)###
######  ###
######  visium <- AddMetaData(visium, metadata = norm_weights)###
######  for (col in colnames(weights)) {###
######    visium@meta.data[[col]][is.na(visium@meta.data[[col]])] <- 0###
######  }###
######browser()###
######  ###
######  return(visium)###
######}###
###
#### UBC cluster order
###clusterOrder <- c(7,2,8,3,0,4,1,6,5)
###cat(sprintf("UBC cluster order: [ %s ]\n",
###    paste(clusterOrder, collapse = ", ")))
###
###mdataCol <- "SCT_snn_res.0.5"
###cat("****************************************************\n")
###cat(sprintf("Deconvoluting against %s\n", mdataCol))
###cat("****************************************************\n")
###for (i in 3:4) {
###  message(sprintf("\tSlide %i", i))
###  outFile <- sprintf("%s/Visium_%i_%s_doubletMode_%s.rds", 
###    outDir, i, mdataCol, dt)
###
###    if (file.exists(outFile)) {
###        cat("Found deconvolution file. Reading\n")
###        cat(sprintf("File exists: %s\n", outFile))
###        t0 <- Sys.time()
###        tmp <- readRDS(outFile)
###        cat("done in ", Sys.time()-t0, "\n")
###    }
###    else {
###        cat("About to deconvolute\n")
###  t0 <- Sys.time()
###  tmp <- deconv_hb_visium_RCTD(
###    visium = visium[[i]],
###    ref = srat, 
###    ref_cell_type = mdataCol,
###    max_cores = 16
###  )
###    cat("done in ", Sys.time()-t0, "\n")
###
###    cat("Saving RDS\n")
###    t0 <- Sys.time()
###    saveRDS(tmp, outFile)
###    cat("done in ", Sys.time()-t0, "\n")
###}
##cat("Plotting\n")
## .plot <- SpatialPlot(
##    tmp, 
##    features = clusterOrder, #unique(as.character(srat[[mdataCol]][,1])), 
##      #  alpha = 0.3, 
##        keep.scale = "feature",
##        min.cutoff = 'q1',
##        max.cutoff = 'q99',
##        label.size=5 
##)
##
##outFile <- sprintf(
##      "%s/Visium%i_%s_%s.png", 
##      outDir, i, mdataCol, dt)
##  ggsave(file=outFile, .plot, 
##    width = 22, height = 20, dpi = 300
##)
##
##  gc()
##}
}, error = function(e) {
  cat("Error: ", e$message, "\n")
  print(e)
}, finally = {
  sink()
})
