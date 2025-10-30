library(Seurat)
rm(list=ls())

inFile <- "/home/rstudio/isilon/public/HumanMouseUBC/data/UBC.Harmony.RDS"

cat("reading file\n")
t0 <- Sys.time()
srat <- readRDS(inFile)
cat("done in ", Sys.time()-t0, "\n")

mdata <- srat[[]]

