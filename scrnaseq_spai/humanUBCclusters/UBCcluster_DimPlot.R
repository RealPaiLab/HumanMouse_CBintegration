
rm(list=ls())
library(Seurat)
library(ggplot2)
library(reshape2)
library(ggpubr)
#srat_file <- "/home/rstudio/isilon/public/HumanMouseUBC/data/UBC.Harmony.RDS"
srat_file <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/DimPlots"

DEdir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/diffExpr"
useDEfile <- sprintf("%s/250424/de_genes.tsv", DEdir)

TFfile <- "/home/rstudio/isilon/src/gene-regulation/Lambert2018_humanTFs/TF_names_v_1.01.txt"

pySCENICFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/20250319/human_ubcs.rss.csv"

colourFile <- "UBCcolours.txt"
clrs <- read.delim(colourFile,header=TRUE,sep="\t")

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outDir,dt)

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

cat("Reading human UBC cluster file\n")
t0 <- Sys.time()
srat <- readRDS(srat_file)
print(Sys.time() - t0)

browser()

xpr <- srat[["RNA"]]$counts

tf <- read.delim(TFfile,header=FALSE,sep="\t")
rss <- read.delim(pySCENICFile,header=TRUE,sep=",",as.is=TRUE)
topregulons <- c()
for (i in 1:nrow(rss)) {
  tmp <- colnames(rss)[which(rss[i,] > 0.3)]
  cat(sprintf("%s: %i\n",rss[i,1],length(tmp)))
  topregulons <- c(topregulons,tmp )
}
topregulons <- setdiff(unique(topregulons),"X")

plotGeneByCluster <- function(geneName) {
  g <- xpr[which(rownames(xpr)==geneName),]
  g <- log2(g+1)
  g <- melt(g)
  x <- cbind(rownames(srat[[]]), 
          as.character(srat[[]]$SCT_snn_res.0.5)
  )
  colnames(x) <- c("Cell","Cluster")
  g$Cell <- rownames(g)
  colnames(g)[1] <- "g"
  y <- merge(x=x,y=g,by="Cell")
  y$Cluster <-  factor(y$Cluster, levels=clrs$Cluster)

  y <- y[-which(y$Cluster %in% c(6)),]

  if (geneName %in% tf[,1]){
    sp <- "(TF)"
  } else {
    sp <- ""
  }

  if (geneName %in% topregulons){
    sp <- paste(sp,"(topReg)")
  }

  p <- ggplot(y,aes(x=Cluster, y=g)) + geom_boxplot(aes(fill=Cluster))
  p <- p + scale_fill_manual(values=clrs$Colour)
  p <- p + ggtitle(sprintf("%s %s", geneName, sp)) +
    xlab("UBC cluster") + ylab("log2(counts+1)")
  p <- p + theme_minimal(base_size = 20)
  p <- p + theme(legend.position = "none")
  #ggsave(p,file=sprintf("%s/boxplot_%s.pdf",outDir,geneName), 
  #    width = 8, height = 6, dpi = 300)

  return(p)
}
#srat$
p <- DimPlot(srat, reduction = "umap", 
    group.by = "SCT_snn_res.0.5", 
    label = TRUE, label.size = 10, pt.size = 1.2,
    label.box=TRUE, repel = TRUE,
    cols=clrs$Colour[order(clrs$Cluster)]) +
  ggtitle(sprintf("Human UBC clusters, %i cells",ncol(srat))) 
p <- p + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave(filename = sprintf("%s/UBCclusters_SCT_snn_res.0.5.pdf", outDir), 
    plot = p, width = 8, height = 6, dpi = 300)

p <- FeaturePlot(srat, features=c("EYS","ONECUT2","MKI67","EOMES", "OTX2"))
ggsave(filename = sprintf("%s/UBCclusters_FeaturePlot.pdf", outDir), 
    plot = p, width = 8, height = 6, dpi = 300)


dge <- read.delim(useDEfile,header=TRUE,sep="\t")
dge2 <- subset(dge, avg_log2FC > 1.2 & p_val_adj < 0.05)


top_genes <- c()
top_gene_regulons <- c()
for (cur in paste("UBC_",clrs$Cluster,sep="")) {
  x <- subset(dge2, ubc_subcluster == cur)
  x <- x[order(x$avg_log2FC,decreasing=TRUE),]
  #if (length(x$gene %in% topregulons)>0){
    #top_gene_regulons <- c(top_gene_regulons, x[which(x$gene %in% topregulons),])
    #x <- x[which(x$gene %in% tf[,1]),]
  #  cat(sprintf("%s: Picking top regulons\n", cur))
  #} else {
    x <- x[1:min(5,nrow(x)),]
  #}
  top_genes <- c(top_genes, x$gene)
  cat(sprintf("%s: %i genes\n", cur, nrow(x)))
}
top_genes <- c(unique(top_genes), c("OTX2","EOMES","BRCA1"))
cat(sprintf("Total top genes = %i\n", length(top_genes)))

##geneList <- c(
##  "MKI67","TOP2A",# 7
##  "ITPR1","SLC6A1", # 3
##  "GREM2", # 8
##  "EYS","IQCJ-SCHIP1", "SEZ6L",#8-3-2
##  "LAMA2","RBFOX1",
##  "OTX2","EOMES","LRP1B") # 0 and 5
##  #"HBG2", "JUND", "PRLR","TACR1","ANKRD18CP","CXCL12" - not distinctive
##  #"SOX11", - "MGAT4C", too many clusters
##  # "PAX2",
##  # "SEMA3A", too many clusters

plotlist <- list()
for (g in top_genes) {
  cat("Plotting gene",g,"\n")
  plotlist[[g]] <- plotGeneByCluster(g)
}

srat$UBC_subcluster <- paste("UBC", srat$SCT_snn_res.0.5,sep="")
lvl <- paste("UBC",clrs$Cluster,sep="")
srat$UBC_subcluster <- factor(srat$UBC_subcluster, levels=lvl)
p <- DotPlot(srat, features = unique(top_genes), group.by = "UBC_subcluster", 
    cols = c("lightgrey", "blue"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = sprintf("%s/UBCclusters_DotPlot.pdf", outDir), 
    plot = p, width = 12, height = 6, dpi = 300)

p <- ggarrange(plotlist = plotlist[1:20])
ggsave(filename = sprintf("%s/UBCclusters_boxplots1.pdf", outDir), 
    plot = p, width = 18, height = 12)

p <- ggarrange(plotlist = plotlist[21:length(top_genes)])
ggsave(filename = sprintf("%s/UBCclusters_boxplots2.pdf", outDir), 
    plot = p, width = 18, height = 12)
