# examine splits of Leo's UBC clusters
library(Seurat)
library(ggpubr)
suppressMessages(library(ggplot2))

setFile <- "/home/rstudio/isilon/private/llau/results/integrated/20240618/rl_cca.qs"
normFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/Integrated_RL_Hsonly/rl_cca_HumanUBCOnly_NormSCTdone_241029.qs"

dt <- format(Sys.Date(),"%y%m%d")

if (!file.exists(normFile)) {
    srat <- qs::qread(setFile)
    browser()
    srat <- subset(srat, species == "human" & broad_cell_type == "UBC")
    srat <- NormalizeData(srat)
    srat <- SCTransform(srat, 
        vars.to.regress="CC.Difference") 
} else {
    srat <- qs::qread(normFile)
}

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters"
outDir <- sprintf("%s/%s",outDir,dt)
if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)

resSet <- c(0.3,0.4,0.6,0.8)
mdata <- srat@meta.data
for (currRes in resSet) {
    message(sprintf("Resolution = %s",currRes))
    cellcats <- mdata[,currRes]

    plist <- list()
    for (k in unique(cellcats)) {
        message(sprintf("\t%s",k))
        srat$cur_cluster <- mdata[,currRes] == k
        plist[[k]] <- DimPlot(srat, reduction="umap", 
            cols=c("grey90","#aa0000"),
            group.by="cur_cluster",
            order=TRUE,
            label=FALSE) + ggtitle(sprintf("%s: Cluster %s",currRes,k)) + 
            xlab("UMAP 1") + ylab("UMAP 2") + 
            theme(plot.title=element_text(size=9)) + NoLegend()
        #ggsave(p,file=sprintf("%s/test.png",outDir))
    }
    p <- ggarrange(plotlist=plist)
    nr <- ceiling(length(unique(cellcats))/3)
    ggsave(p,
        file=sprintf("%s/DimPlot_%s_%s.png",
            outDir,currRes,dt),
        width=8,height=nr*3,unit="in"
        )
}