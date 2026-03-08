
rm(list=ls())
library(Seurat)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

srat_file <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/humanDEG"

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE);  
}

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE);  
}

cat("Reading human UBC cluster file\n")
t0 <- Sys.time()
srat <- readRDS(srat_file)
print(Sys.time() - t0)

# remove cells with age "9 PCW" and "10 PCW"
srat <- subset(srat, subset = age != "9 PCW" & age != "10 PCW")
cat(sprintf("After removing 9 and 10 PCW, %i cells remain\n", ncol(srat)))
srat_mega <- srat

srat <- subset(srat_mega, PAX6 > 0 & EOMES > 0)
cat(sprintf("After subsetting to PAX6+ EOMES+ cells, %i cells remain\n", ncol(srat)))
Idents(srat) <- srat$SCT_snn_res.0.5

p <- FeaturePlot(srat, features = c("MKI67","CST3", "LRP1B", "CNTNAP2"), order = TRUE, pt.size = 1.2)
ggsave(filename = sprintf("%s/Featureplot.pdf", outDir), plot = p, width = 8, height = 6, dpi = 300)