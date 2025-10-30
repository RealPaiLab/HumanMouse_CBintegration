rm(list=ls())

library(tidyverse)
library(patchwork)
library(Seurat)
library(pals)
library(readxl)


###source("FromIan/utils.R")
###source("FromIan/utilities/plotting.R")
###source("FromIan/utilities/cluster_barplot.R")
###source("FromIan/utilities/plot_venn_diagrams.R")
###source("SP_utils.R")
#source("utilities/propeller_helpers.R")

oldRoot <- "/.mounts/labs/pailab"
newRoot <- "/home/rstudio/isilon" # docker mount point

evoFile <- "/home/rstudio/isilon/src/evolution/Florio_Huttner_2015/aaa1975_supportingfile_suppl.excel_seq2_v1_tables2.xlsx"

out_dir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBC1_marker_analysis"

dt <- format(Sys.Date(),"%y%m%d")

huttner57 <- read_excel(evoFile,sheet=4) 
huttner57 <- na.omit(huttner57$external_gene_id)

# hendrikse dataset
# /.mounts/labs/pailab/private/icheong/CBL_scRNAseq/data/human/Aldinger/glutamatergic_dev_Liam.RDS

datasetFile <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/scrnaseq_Leo/integrations/dataset_list.csv"
dat <- read.delim(datasetFile,sep=",")

groupByList <- list(
    Aldinger_RL_human="common_cell_name",
    Sepp_RL_human="common_cell_name",
    Luo_RL_human="common_cell_name"
)


#setName <- "Luo_RL_human"
#setName <- "Aldinger_RL_human"
setName <- "Integrated_RL"
groupBy <- "common_cell_name"

if (setName %in% "Sepp_RL_human") {
    setFile <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20241017/SEPP_RL_human_norm.qs"
} else if (setName %in% "Integrated_RL") {
    setFile <- "/home/rstudio/isilon/private/llau/results/integrated/20240618/rl_cca.qs"
} else {
  setFile <- sub(pattern = oldRoot, replacement = newRoot, 
    x = dat$dataset_location_qs[which(dat[,1] %in% setName)])
}

setDir <- sprintf("%s/%s",out_dir,setName)
if (!file.exists(setDir)) dir.create(setDir)

browser()

message(sprintf("Loading %s", setName))
t0 <- Sys.time()
srat <- qs::qread(setFile)
t1 <- Sys.time()
print(t1-t0)

if (!setName %in% c("Sepp_RL_human")) {
  srat <- UpdateSeuratObject(srat)
  if (setName %in% "Integrated_RL"){
      message("\tSubsetting for human cells")
      srat <- subset(rl_srat, subset = species == "human")
  }
    t0 <- Sys.time()
    srat <- NormalizeData(srat)
    t1 <- Sys.time()
    print(t1-t0)
    t0 <- Sys.time()
    srat <- SCTransform(srat, 
    vars.to.regress="CC.Difference") 
    t1 <- Sys.time()
    print(t1-t0)
} 
srat <- RunPCA(srat, verbose = FALSE)
srat <- RunUMAP(srat, dims = 1:30, verbose = FALSE)
srat <- FindNeighbors(srat, dims = 1:30, verbose = FALSE)
srat <- FindClusters(srat, verbose = FALSE)

p <- DimPlot(srat, label = TRUE, group.by=groupBy)
ggsave(p,file=sprintf("%s/DimPlot_%s.png",setDir,dt))

huttner57 <- intersect(huttner57, Features(srat))
cat(sprintf("%i Huttner genes in feature set",length(huttner57)))
geneSet <- c("OTX2","EOMES","ARHGAP11B", "C2orf48","SOX2","MKI67",
"SOX4","SOX11")

browser()

numPages <- ceiling(length(geneSet)/8)
for (k in 1:numPages){
  sidx <- (8*(k-1))+1
  eidx <- min(8*k, length(geneSet))
  message(sprintf("%i: %i - %i",k, sidx,eidx))
  p <- FeaturePlot(srat, geneSet[sidx:eidx],order=TRUE, ncol=4)
  ggsave(p, file=sprintf("%s/FeaturePlot_page%i_%s.pdf",setDir,k,dt),
    width=11, height=6,unit="in")
  
  p <- DotPlot(srat, geneSet[sidx:eidx], assay="SCT",group.by=groupBy)
  ggsave(p, 
    file=sprintf("%s/DotPlot_page%i_%s.png",setDir,k,dt),
    width=11,height=8,unit="in")
}

p <- VlnPlot(srat, features=geneSet, group.by=groupBy)
ggsave(p, file=sprintf("%s/VlnPlot_%s.png",setDir,dt))