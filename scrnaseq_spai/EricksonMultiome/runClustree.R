# run clustree on various tumour resolutions
# this is run on the HPC

require(Seurat)
library(clustree)
mb <- "/home/rstudio/isilon/private/projects/MB_multiome/output/QC/20250206/fastmnn/mb_fastmnn.qs"
#mb <- "/.mounts/labs/pailab/private/projects/MB_multiome/output/QC/20250206/fastmnn/mb_fastmnn.qs"

outDir <- "/.mounts/labs/pailab/private/projects/MB_multiome/output/clustering"
clusterFile <- sprintf("%s/250502/clustering_resolutions.txt", outDir)
outDir <- dirname(clusterFile)
#dt <- format(Sys.Date(),"%y%m%d")

cat("Loading Erickson MB data - takes 13 min on Pai lab server\n")
t0 <- Sys.time()
srat <- qs::qread(mb)
print(Sys.time()-t0)

cat("Now adding cluster assignments to main Seurat\n")
cl <- read.table(clusterFile, header=TRUE,sep="\t")
cols <- colnames(cl)[grep("^snn_res", colnames(cl),ignore.case=FALSE)]
cols <- setdiff(cols, c("snn_res.0.5","snn_res.0.7","snn_res.0.9"))

md <- srat[[]]
if (all.equal(rownames(md),rownames(cl))!=TRUE){
    cat("Row names of metadata and clustering do not match. Fix this. \n")
    browser()
}

for (i in cols) {
    cat("Adding ",i," to metadata\n")
    srat[[i]] <- cl[[i]]
}

srat[["snn_res.0.8"]] <- cl[["SCT_snn_res.0.8"]]

browser()

# Extract and visualize the clusters
cat("\tGenerating clustree plot\n")
t0 <- Sys.time()
tree <- clustree(srat, prefix = "snn_res.")
print(Sys.time() - t0)
ggsave(tree, filename = "Clustree.png", path = outDir, 
            width = 18, height = 12, units = "in")
    