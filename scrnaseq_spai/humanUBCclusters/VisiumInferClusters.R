# Run RCTD on cleaned UBC cluster set with cleaned marker genes
rm(list=ls())

UBC_Seurat <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"
visiumFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/Aldinger_Visium_2021Data/obj.list.rds"

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/HumanUBC/Aldinger_Visium"
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)

library(spacexr)

markerGeneFile<- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/UBCmarkergenes_251002.txt"


#' deconvolute hindbrain visium against scRNA ref using RCTD
#' @param visium (SeuratObject) visium spatial object
#' @param ref (SeuratObject) scRNA reference
#' @param ref_cell_type (character) Column to be used as reference cell type
#' @param max_cores (numeric) Max number of cores for RCTD
#' @return (SeuratObject) The original visium object with deconvoluted cell ###type weigths
deconv_hb_visium_RCTD <- function(visium, ###
                                  ref, ref_cell_type, ###
                                  max_cores = 8) {###
  # extract information to pass to the RCTD Reference function
  counts <- ref@assays$RNA$counts###
  cluster <- as.factor(ref@meta.data[[ref_cell_type]])
  names(cluster) <- rownames(ref@meta.data)
  nUMI <- ref@meta.data$nCount_RNA
  names(nUMI) <- rownames(ref@meta.data)
  reference <- Reference(counts, cluster, nUMI)
  #
  # set up query with the RCTD function SpatialRNA
  counts <- visium[["Spatial"]]$counts
  coords <- GetTissueCoordinates(visium)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))
  #
  # RCTD
  RCTD <- create.RCTD(query, reference, max_cores = max_cores)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  #
  # assign cell type weights
  weights <- RCTD@results$weights
  norm_weights <- normalize_weights(weights)
  visium <- AddMetaData(visium, metadata = norm_weights)
  for (col in colnames(weights)) {
    visium@meta.data[[col]][is.na(visium@meta.data[[col]])] <- 0
  }
  ##browser()
  return(visium)
}

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = T)
}

logFile <- sprintf("%s/log.txt", outDir)
sink(logFile, split=TRUE)

tryCatch({

cat("Reading UBC marker genes ...\n")
markerGenes <- read.delim(markerGeneFile, 
    header=FALSE, stringsAsFactors=FALSE)[,1]
markerGenes <- unique(markerGenes)
cat(sprintf("Read %i unique marker genes from %s\n", 
    length(markerGenes), markerGeneFile))
cat("removing ribosomal proteins, which are highly expressed, and adding EOMES")
markerGenes <- c(markerGenes[!grepl("^RPS|^RPL", markerGenes)], "EOMES")
cat(sprintf("keeping %i marker genes\n", length(markerGenes)))

cat("Reading UBC clusters ...")
cat(sprintf("Input file:\n%s\n", UBC_Seurat))
t0 <- Sys.time()
srat <- readRDS(UBC_Seurat);
cat("done in ", Sys.time()-t0, "\n")

DefaultAssay(srat) <- "RNA"

cat("exclude ribosomal proteins from Seurat object\n")
gn <- rownames(srat)
tokeep_idx <- !grepl("^RPS|^RPL", gn)
srat <- srat[tokeep_idx,]

#srat2 <- subset(srat, features = markerGenes)

clusterOrder <- c(7,2,8,3,0,4,1,5)

cat(sprintf("Filtered to %i genes\n", nrow(srat)))
cat("Filtering for high confidence UBC clusters")
#srat <- srat3

mdataCol <- "SCT_snn_res.0.5"  # UBC clusters

# UBC cluster order

cat(sprintf("UBC cluster order: [ %s ]\n",
    paste(clusterOrder, collapse = ", ")))

cat("Reading Visium data\n")
visium <- readRDS(visiumFile)
    lbls <- c("sample1_slice1", "sample1_slice2", "sample2_slice1", "sample2_slice2")

cat("****************************************************\n")
cat(sprintf("Deconvoluting against %s\n", mdataCol))
cat("****************************************************\n")
for (i in 3) {
  message(sprintf("\tSlide %i", i))
  outFile <- sprintf("%s/Visium_%i_%s_doubletMode_%s.rds", 
    outDir, i, mdataCol, dt)

    if (file.exists(outFile)) {
        cat("Found deconvolution file. Reading\n")
        cat(sprintf("File exists: %s\n", outFile))
        t0 <- Sys.time()
        tmp <- readRDS(outFile)
        cat("done in ", Sys.time()-t0, "\n")
    }
    else {
        cat("About to deconvolute\n")
  t0 <- Sys.time()
  tmp <- deconv_hb_visium_RCTD(
    visium = visium[[i]],
    ref = srat, 
    ref_cell_type = mdataCol,
    max_cores = 16
  )
    cat("done in ", Sys.time()-t0, "\n")

    cat("Saving RDS\n")
    t0 <- Sys.time()
    saveRDS(tmp, outFile)
    cat("done in ", Sys.time()-t0, "\n")
}
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
}

}, error = function(e) {
  cat("Error: ", e$message, "\n")
  print(e)
}, finally = {
  sink()
})