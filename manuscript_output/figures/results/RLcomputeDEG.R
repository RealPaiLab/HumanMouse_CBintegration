# compute DEG for rhombic lip cells in development
rm(list=ls())

library(argparse)
library(tidyverse)
library(patchwork)
library(Seurat)
library(SingleR)

ref_rds <- "/home/rstudio/isilon/private/llau/results/integrated/20240524/25_pc_without_luo/25_pc_rl.qs"
out_dir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/RLDEG"

if (!file.exists(out_dir)) {
    dir.create(out_dir, recursive = FALSE)
}

dt <- format(Sys.Date(), "%y%m%d")
out_dir <- sprintf("%s/%s", out_dir, dt)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = FALSE)
}


logFile <- sprintf("%s/RL_DEG.txt", out_dir)
sink(logFile, split=TRUE)
tryCatch({

# reference
cat("Reading reference...")
t0 <- Sys.time()
ref_srat <- qs::qread(ref_rds)
print(Sys.time() - t0)
cat("done\n")

ref_srat <- subset(ref_srat, species == "human")
cat(sprintf("Number of reference cells after subsetting for human: %d\n", ncol(ref_srat)))

browser()

srat <- ref_srat

cat("Running SCTransform\n")
options(future.globals.maxSize= 50*1024^3)
srat$percent.mt <- PercentageFeatureSet(srat, pattern = "^MT-")
t0 <- Sys.time()
srat <- SCTransform(
    srat, 
    vars.to.regress = c("CC.Difference","percent.mt","dataset_name"),
    return.only.var.genes = FALSE,
    verbose = FALSE
)
print(Sys.time() - t0)

srat <- PrepSCTFindMarkers(srat)
srat$seurat_clusters <- factor(srat$common_cell_name)
source("clusterUtils.R")
de <- run_all_de(
    srat,
    clusters="RL",
    grouping.var = "dataset_name",
    logfc.threshold = 0.2,
    min.pct=0.2
)

write.table(de, file=sprintf("%s/RL_DEG.txt", out_dir), sep="\t", quote=FALSE, row.names=FALSE)

p <- plotVolcano_usingDEG(de, showTopGenes=50, titlePfx="RL")
ggsave(p, file=sprintf("%s/RL_DEG_volcano.pdf", out_dir), width=8, height=6)

}, error=function(ex){
    print(ex)
}, finally={
    
    sink()
    message("\n***SESSION INFO***\n")
    print(sessionInfo())

})




