# use featureplot to look at cluster integrity of ubc clusters one by one
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(Seurat)

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBC1_marker_analysis/Integrated_RL_Hsonly"
normFile <- sprintf("%s/rl_cca_HumanOnly_NormSCTdone_241029.qs",
    outDir)
dt <- format(Sys.Date(),"%y%m%d")

srat <- qs::qread(normFile)
srat <- RunUMAP(srat, dims = 1:25, verbose = FALSE)    
human_cells <- length(Cells(srat))
p <- DimPlot(srat, reduction = "umap", group.by = "broad_w_ubc_subtypes", 
    label = TRUE, label.size = 4, raster = FALSE, repel = TRUE) 
ggsave(p,file=sprintf("%s/DimPlot_%s.png",outDir,dt),
    width=8, height=6, unit="in")

plist <- list()
cellcat <- c(
    "RL","GCP","GC",
    sprintf("UBC_%i",0:5),
    "microglia","endothelial","oligodendrocyte/OPC"
)
for (k in cellcat) {
    srat$cur_cluster <- mdata$broad_w_ubc_subtypes == k
    plist[[k]] <- DimPlot(srat, reduction="umap", 
        cols=c("grey90","#aa0000"),
        group.by="cur_cluster",
        order=TRUE,
        label=FALSE) + ggtitle(k) + 
        xlab("UMAP 1") + ylab("UMAP 2") + 
        theme(plot.title=element_text(size=9)) + NoLegend()
    #ggsave(p,file=sprintf("%s/test.png",outDir))
}
p <- ggarrange(plotlist=plist)
ggsave(p,
    file=sprintf("%s/Allclusters_DimPlot_Single_%s.pdf",
        outDir,dt),
    width=8,height=8,unit="in"
)

p <- ggarrange(plotlist=plist[grep("UBC",names(plist))])
ggsave(p,
    file=sprintf("%s/UBCclusters_DimPlot_Single_%s.pdf",
        outDir,dt),
    width=8,height=5,unit="in"
)

