# plot expression of gene set in UBC 1 & 2

# load Aldinger & Sepp
# batch correct
# subset for cells in UBC1 and 2
# plot their xpr

rm(list=ls())

library(tidyverse)
library(patchwork)
library(Seurat)
library(pals)
library(readxl)
library(ggpubr)
oldRoot <- "/.mounts/labs/pailab"
newRoot <- "/home/rstudio/isilon" # docker mount point

out_dir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBC1_marker_analysis"
dt <- format(Sys.Date(),"%y%m%d")
ubcID <- sprintf("%s/humanUBC_CellID.txt",out_dir)
ubcID <- read.delim(ubcID,sep="\t",h=T,as.is=T)
ubcID <- ubcID[,c("subclusters","cell_id")]
ubcID[,2] <- as.character(ubcID[,2])

datasetFile <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/scrnaseq_Leo/integrations/dataset_list.csv"
dat <- read.delim(datasetFile,sep=",")

groupByList <- list(
    Aldinger_RL_human="common_cell_name",
    Sepp_RL_human="common_cell_name",
    Luo_RL_human="common_cell_name"
)

#setName <- "Aldinger_RL_human"
setName <- "Aldinger_RL_human"
groupBy <- "common_cell_name"

if (setName %in% "Sepp_RL_human") {
    setFile <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20241017/Sepp_RL_human_norm.qs"
} else {
setFile <- sub(pattern = oldRoot, replacement = newRoot, 
  x = dat$dataset_location_qs[which(dat[,1] %in% setName)])
}

setDir <- sprintf("%s/%s",out_dir,setName)
if (!file.exists(setDir)) dir.create(setDir)

message(sprintf("Loading %s", setName))
t0 <- Sys.time()
srat <- qs::qread(setFile)
t1 <- Sys.time()
print(t1-t0)

DefaultAssay(srat) <- "RNA"
if (!setName %in% c("Sepp_RL_human")) {
    srat <- UpdateSeuratObject(srat)
    message("Normalizing")
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

geneSet <- c("FOXP2","EOMES","OTX2","ARHGAP11B","C2orf48","CNTNAP2","RBFOX1",
"SOX5","MYT1L", "EBF1","ST18","SOX2","KI67","SOX4","SOX11")

geneSet <- sort(intersect(geneSet,Features(srat)))
browser()

plotlist <- list()
RLlist <- list()
for (geneName in geneSet) {
    print(geneName)
    xpr <- srat[["SCT"]]$counts[geneName,]
    xpr <- as.data.frame(xpr)
    xpr$cell_id <- as.character(rownames(xpr))

    x <- merge(x=xpr,y=ubcID,by="cell_id" )
    x$subclusters[which(!x$subclusters %in% "UBC_1")] <- "other_UBC"
    x$subclusters <- factor(x$subclusters, levels=c("UBC_1","other_UBC"))

    pval <- t.test(x$xpr[which(x$subclusters %in% "UBC_1")],
        x$xpr[which(!x$subclusters %in% "UBC_1")])$p.value
    
    plotlist[[geneName]] <- ggplot(x,aes(x=subclusters, y=xpr)) + 
        geom_violin(aes(fill=subclusters)) + 
        geom_point(size=2) +
        ggtitle(sprintf("%s (p < %1.2e)", geneName, pval))
}

megap <- ggarrange(plotlist=plotlist)
ggsave(megap, 
    file=sprintf("%s/UBC_boxplots.png",setDir),
    width=11,height=8,unit="in"
)

##p <- DotPlot(srat, geneSet, group.by="subclusters")
##ggsave(p,
##    file=sprintf("%s/UBC_dotplots.png",setDir)
##)


