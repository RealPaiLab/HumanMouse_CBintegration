#rm(list=ls())
library(Seurat)
# run DEG analysis
#mb <- "/home/rstudio/isilon/private/projects/MB_multiome/output/QC/20250206/fastmnn/mb_fastmnn.qs"
mb <- "/.mounts/labs/pailab/private/projects/MB_multiome/output/QC/20250206/fastmnn/mb_fastmnn.qs"

outRoot <- "/.mounts/labs/pailab/private/projects/MB_multiome/output"
clusterFile <- sprintf("%s/clustering/250502/clustering_assignments.txt", outRoot)
phenoFile <- "/.mounts/labs/pailab/private/projects/MB_multiome/input/AndersErickson/20241212_multiome_metadata.tsv"

devProjections <- "/.mounts/labs/pailab/private/projects/HumanMouseUBC/SingleR_UBCclusters_Erickson/250506/Erickson_singleRPred_250506.RData"

outDir <- sprintf("%s/DEGanalysis_byDevType",outRoot)
dt <- format(Sys.Date(),"%y%m%d")
outDir <- file.path(outDir,dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE)
}

logFile <- sprintf("%s/runDE.log", outDir)
sink(logFile,split=TRUE)
tryCatch({

 # read pheno file
cat("Reading pheno file\n")
pheno <- read.table(phenoFile, header=TRUE, sep="\t")

 cat("Loading Erickson MB data - takes 13 min on Pai lab server\n")
    t0 <- Sys.time()
  #  srat <- qs::qread(mb)
    print(Sys.time()-t0)

browser()
cat("Integrating sex info\n")
mdata <- srat[[]]
blah <- pheno[,c("unique_id","sex")]
midx <- match(mdata$sample_id, blah$unique_id)
if (all.equal(blah$unique_id[midx],mdata$sample_id)!=TRUE){
    cat("Row names of metadata and pheno do not match. Fix this. \n")
    browser()
}
srat$sex <- blah$sex[midx]

cat("Now adding cluster assignments to main Seurat\n")
cl <- read.table(clusterFile, header=TRUE,sep="\t")

md <- srat[[]]
if (all.equal(rownames(md),rownames(cl))!=TRUE){
    cat("Row names of metadata and clustering do not match. Fix this. \n")
    browser()
}
srat[["snn_res.1"]] <- factor(cl[,"snn_res.1"])

cat("Adding developmental projections\n")
load(devProjections)

browser()

mdata <- srat[[]]
mdata2 <- mdata[,c("subtype","seurat_clusters")]
mdata2$CellID <- rownames(mdata2)

singleRPred$CellID <- rownames(singleRPred)
merged <- merge(mdata2, singleRPred, by = "CellID")

cat("Add SingleR labels to srat\n")
midx <- match(rownames(mdata), merged$CellID)
if (all.equal(merged$CellID[midx],rownames(mdata))!= TRUE) {
    stop("CellID mismatch")
}

srat$RLdev_SingleR_pruned.labels <- merged$pruned.labels[midx]
srat$RLdev_SingleR_labels <- merged$labels[midx]


cat("Filtering for genes expressed in >1%% of cells\n")
# If a count for a gene in a cell is greater than 0, set as TRUE (= 1)
counts <- GetAssayData(object = srat, slot = "counts")
nonzero <- counts > 0
# If 1% or more cells are TRUE, keep the gene. Each TRUE value = 1. Taking the sum of all the cells for that gene
cell_count <- Matrix::rowSums(nonzero) 
keep_genes <- which(cell_count >= (0.01*length(Cells(srat))))

cat(sprintf("\t\t%i total, %i in 1%%\n", length(rownames(srat)),
length(keep_genes)))
# Get the genes names that we are keeping
keep_genes <- rownames(srat)[keep_genes]

rm(cell_count)
rm(nonzero)
rm(counts)
gc()

#full_srat <- srat
t0 <- Sys.time()
srat <- subset(srat, features = keep_genes)
print(Sys.time()-t0)

cat("\t\tPrepSCTFindMarkers\n")
t0 <- Sys.time()
srat <- PrepSCTFindMarkers(object = srat)
print(Sys.time()-t0)

cat("Just finished PrepSCTFindMarkers\n")
browser()
pattern <- "RLdev_SingleR_labels"

# Setting the Idents to be the cluster name because FindMarkers() uses the Ident
Idents(srat) <- pattern

# remove samples with sex NA
cat(sprintf("\t\t%i cells before\n", length(Cells(srat))))
cat("\t\tRemoving samples with NA sex\n")
srat <- subset(srat, subset = !is.na(srat$sex))
cat(sprintf("\t\t%i cells remaining\n", length(Cells(srat))))
browser()

cat("\t\tRunning DEG analysis\n")
marker_list <- list()
# Find markers for each cluster
for (cluster in levels(Idents(srat))) {
    print(cluster)
     markers <- FindMarkers(
        object = srat,
        ident.1 = cluster,
        assay = "SCT",
        logfc.threshold = 0.7,
        min.pct = 0.1,
        test.use = "MAST",
        latent.vars = c("sex")
      ) 

    cat("writing markers\n")
    write.table(markers,
        file=sprintf("%s/%s_DEG_sexCovar.csv",outDir,cluster),
        sep=",",quote=FALSE)
    
# Append markers to list
    marker_list[[cluster]] <- markers
}
browser()

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
}, error = function(e) {
    print(e)
}, finally = {
    cat("Closing log.\n")
	print(sessionInfo())
    sink(NULL)
})
