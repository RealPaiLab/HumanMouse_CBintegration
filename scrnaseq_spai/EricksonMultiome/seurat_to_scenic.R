# ==============================================================================
# Save Seurat raw counts and metadata as CSV for use in SCENIC.
# ==============================================================================

library(argparse)
library(tidyverse)
library(Seurat)

inFile <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/tumour/MB_multiome/20250206/fastmnn/mb_fastmnn.qs"
outDir <- "/home/rstudio/isilon/private/projects/MB_multiome/output/PySCENIC"

baseF <- basename(inFile)
countFile <- sprintf("%s/%s_counts.txt", outDir, baseF)
metadataFile <- sprintf("%s/%s_metadata.txt", outDir, baseF)
umapFile <- sprintf("%s/%s_umap.txt", outDir, baseF)

dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}


cat("Reading tumour file\n")
t0 <- Sys.time()
srat <- qs::qread(inFile)
print(Sys.time() - t0)

counts <- GetAssayData(srat[["RNA"]], slot = "counts")
cat("writing counts\n")
browser()
ct <- as.data.frame(counts)
write.csv(ct, file = countFile)

message(sprintf("Saving metadata to %s", metadataFile))
write.csv(srat[[]], file = metadataFile)

cat("saving UMAP\n")
umap_embeds <- srat[["umap"]]@cell.embeddings %>% as.data.frame()
write.csv(umap_embeds, file = umapFile)

message("\nSESSION INFO\n")
print(sessionInfo())

