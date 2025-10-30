# compute UBC1 module score

rm(list=ls())

library(tidyverse)
library(patchwork)
library(Seurat)
library(pals)
library(ggpubr)

out_dir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBC1_marker_analysis"
dt <- format(Sys.Date(),"%y%m%d")


oldRoot <- "/.mounts/labs/pailab"
newRoot <- "/home/rstudio/isilon" # docker mount point
datasetFile <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/scrnaseq_Leo/integrations/dataset_list.csv"
dat <- read.delim(datasetFile,sep=",")

MBfile <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20230510/mb_mnn.rds"
th <- theme(text=element_text(size=20))

setName <- "Integrated_RL"
groupBy <- "broad_w_ubc_subtypes"
setDir <- sprintf("%s/%s",out_dir,setName)
setFile <- "/home/rstudio/isilon/private/llau/results/integrated/20240618/rl_cca.qs"

# ----------------------------------
message("Generate UBC1 signature")
# load UBC1 genes
DEGresults <- "/home/rstudio/isilon/private/llau/results/integrated/20240715/all_tested_genes.csv"

dat <- read.delim(DEGresults,sep=",",h=T,as.is=T)
dat <- subset(dat,cluster %in% "UBC_1")
dat <- subset(dat, 
    Aldinger_full_cerebellum_human_p_val_adj < 0.05 & 
    Sepp_full_cerebellum_human_p_val_adj < 0.05
)

dat <- dat[,c(1,2,6,7,14)]
dat$avgLogFC <- (dat[,2]+dat[,4])/2
dat$Ald <- dat$Aldinger_full_cerebellum_human_avg_log2FC
dat$Sepp <- dat$Sepp_full_cerebellum_human_avg_log2FC
dat <- dat[order(dat$avgLogFC,decreasing=TRUE),]
cat(sprintf("UBC1 vs other UBC\n"))
cat(sprintf("%i genes with Q < 0.05 in Aldinger & Sepp\n",nrow(dat)))

ubc1 <- dat$feature[1:20]
cat(sprintf("Signature of top 20 genes = { %s }\n", 
    paste(ubc1,collapse=",")))

# ----------------------------------
message("Project onto developmental dataset")
message(sprintf("Loading %s", setName))
t0 <- Sys.time()
srat <- qs::qread(setFile)
srat_orig <- srat
t1 <- Sys.time()
print(t1-t0)

srat <- UpdateSeuratObject(srat)
srat <- subset(srat, subset = species == "human")
srat <- NormalizeData(srat, assay="RNA")
srat <- RunUMAP(srat, dims = 1:25, verbose = FALSE)    
human_cells <- length(Cells(srat))

p <- DimPlot(srat, reduction = "umap", group.by = "broad_w_ubc_subtypes", 
    label = TRUE, label.size = 4, raster = FALSE, repel = TRUE) +
    ggtitle(paste0("Rhombic Lip Cell Types ", "(n = ", human_cells, ")")) 
ggsave(p,file=sprintf("%s/DimPlot_%s.png",setDir,dt),
    width=8, height=6, unit="in")

srat <- AddModuleScore(srat,
    features = list(ubc1),
    name="UBC1_sig")

p <- FeaturePlot(srat, 
    features="UBC1_sig1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = 
        rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(p,file=sprintf("%s/UBC1score_%s.png",setDir,dt),
    width=8, height=6, unit="in")

xpr <- srat@meta.data[,c("broad_w_ubc_subtypes","UBC1_sig1")]
p <- ggplot(xpr, 
    aes(x=broad_w_ubc_subtypes, y=UBC1_sig1))+
    geom_violin(aes(fill=broad_w_ubc_subtypes))+
    geom_boxplot(aes(fill=broad_w_ubc_subtypes),width=0.1) + 
    xlab("Major cell type") + ylab("UBC1 signature") # +
    #geom_point(aes(fill=subtype))
ggsave(p,file=sprintf("%s/UBC1score_Violin_by_CellType_%s.png",setDir,dt),
    width=11, height=5)

# ----------------------------------
# load Vladoiu tumour set.
message("Project onto MB dataset")
mb <- readRDS(MBfile)
p <- DimPlot(mb, group.by="seurat_clusters")
ggsave(p,file=sprintf("%s/Vladoiu_MB2019/DimPlot_bycluster_%s.png",
    out_dir,dt))
p <- DimPlot(mb, group.by="subtype")
ggsave(p,file=sprintf("%s/Vladoiu_MB2019/DimPlot_bysubtype_%s.png",
    out_dir,dt))

mb <- AddModuleScore(mb,
    features=list(ubc1),
    name="UBC1_sig")
p <- FeaturePlot(mb,
    features="UBC1_sig1",label=TRUE, repel=TRUE) +
    scale_colour_gradientn(colours = 
        rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(p,file=sprintf("%s/Vladoiu_MB2019/UBC1_sig_%s.png",out_dir,dt))

p1 <- VlnPlot(mb, features="UBC1_sig1", group.by="subtype")
p2 <- VlnPlot(mb, features="UBC1_sig1", group.by="seurat_clusters")
plist <- ggarrange(plotlist=list(p1,p2),nrow=2)
ggsave(plist, 
    file=sprintf("%s/Vladoiu_MB2019/UBC1_sig_violin_%s.png",
        out_dir,dt))

xpr <- mb@meta.data[,c("seurat_clusters","subtype","UBC1_sig1")]

# plot signature levels different tumour clusters
p <- ggplot(xpr, 
    aes(x=subtype, y=UBC1_sig1))+
    geom_violin(aes(fill=subtype))+
    geom_boxplot(aes(fill=subtype),width=0.1) + 
    xlab("Tumour subgroup") + ylab("UBC1 signature") + th # +
    #geom_point(aes(fill=subtype))

ggsave(p,
    file=sprintf("%s/Vladoiu_MB2019/UBC1_sig_violin_bysubtype_%s.png",
    out_dir,dt))

message("Compare signature levels in G4 vs G3 and SHH tumours")
pval <- t.test(xpr[which(xpr$subtype=="G4"),3], 
    xpr[which(xpr$subtype=="G3"),3], alternative="greater")
pval <- t.test(xpr[which(xpr$subtype=="G4"),3], 
    xpr[which(xpr$subtype=="SHH"),3], alternative="greater")


# plot signature levels different tumour clusters
p <- ggplot(xpr, 
    aes(x=seurat_clusters, y=UBC1_sig1))+
    geom_violin(aes(fill=seurat_clusters))+
    geom_boxplot(aes(fill=seurat_clusters),width=0.1)  + th # +
    #geom_point(aes(fill=subtype))

ggsave(p,
    file=sprintf("%s/Vladoiu_MB2019/UBC1_sig_violin_byclusters_%s.png",
    out_dir,dt),width=11,height=5,unit="in")
