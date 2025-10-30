rm(list=ls())

library(Seurat)
library(ggplot2)
library(ggpubr)

inDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/assignClusterIdentity"
dt <- format(Sys.Date(),"%y%m%d")

inFile <- sprintf("%s/cca_RLlineage_only_241107.qs",inDir)
ttl <- "Ald-Sepp CCA > rhombic lip + controls"

outDir <- sprintf("/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/RLlineage_only/%s",dt)

if (!file.exists(outDir)) dir.create(outDir)

geneSet <- c("MKI67","WLS","EOMES","LMX1A","OTX2","RBFOX3","ATOH1","PAX6",
    "RELN","PDGFRA","TNR","TREM2","CLDN5","ITM2A")

UBCmarkers <- c("EOMES","MKI67","CRX","NRL","EYS","CNTNAP2","CALR","ARHGAP11B")

cat("Reading file")
srat <- qs::qread(inFile)
browser()
srat <- PrepSCTFindMarkers(srat)
DefaultAssay(srat) <- "integrated"
srat <- RunPCA(srat, verbose = FALSE)

cat("finding clusters")

for (x in c(15)){
    print(x)
    plist <- list()
    srat <- RunUMAP(srat, dims = 1:x, verbose = FALSE)
    srat <- FindNeighbors(srat, dims = 1:x, verbose = FALSE)
    srat <- FindClusters(srat, verbose = FALSE)

    p <- DimPlot(srat,label=TRUE)
    p <- p + ggtitle(sprintf("%i dims",x)); 
    outFile <- sprintf("%s/%idims_DimPlot_%s.pdf",outDir,x,dt)
    ggsave(p,file=outFile)

    p <- DotPlot(srat, geneSet, assay="SCT") + 
    ggtitle(ttl)
    outFile <- sprintf("%s/%idims_DotPlot_%s.pdf",outDir,x,dt)
    ggsave(p,file=outFile,width=13,height=6,unit="in")
    
}

x <- 15
srat <- RunUMAP(srat, dims = 1:x, verbose = FALSE)
srat <- FindNeighbors(srat, dims = 1:x, verbose = FALSE)
srat <- FindClusters(srat, verbose = FALSE)

p <- DimPlot(srat, label = TRUE)
p <- p + ggtitle(sprintf("%s: cell clusters",ttl))
ggsave(p,file=sprintf("%s/DimPlot_%s.png",outDir,dt))

p <- DimPlot(srat, group.by="sex")
p <- p + ggtitle(sprintf("%s: sex",ttl))
ggsave(p,file=sprintf("%s/DimPlot_sex_%s.png",outDir,dt))

p <- DimPlot(srat, group.by="dataset_name")
p <- p + ggtitle(sprintf("%s: dataset name",ttl))
ggsave(p,file=sprintf("%s/DimPlot_datasetname_%s.png",outDir,dt))


p <- DimPlot(srat, group.by="age")
p <- p + ggtitle(sprintf("%s: age",ttl))
ggsave(p,file=sprintf("%s/DimPlot_age_%s.png",outDir,dt))

p <- DotPlot(srat, geneSet, assay="SCT") + 
    ggtitle(ttl)
ggsave(p,file=sprintf("%s/DotPlot_%s.pdf",outDir,dt),width=13,height=6,unit="in")

srat_RL <- srat; 
# now let's recluster UBC

srat <- subset(srat, seurat_clusters %in% c(0,5,6,11))
x <- 15
srat <- RunPCA(srat,verbose=FALSE)
srat <- RunUMAP(srat, dims = 1:x, verbose = FALSE)
srat <- FindNeighbors(srat, dims = 1:x, verbose = FALSE)

# remove old cluster identities
idx <- grep("snn_",colnames(srat@meta.data))
srat@meta.data[,idx] <- NULL

for (curRes in seq(0.2,0.6,0.2)){
    print(curRes)
    srat <- FindClusters(srat, resolution=curRes, verbose = FALSE)
    mdata <- srat@meta.data
    if (curRes < 1 ) {
        groupBy <- sprintf("integrated_snn_res.%1.1f",curRes)
    } else {
        groupBy <- sprintf("integrated_snn_res.%i",curRes)
    }
    
    p <- DimPlot(srat,group.by=groupBy, label=TRUE)
    p <- p + ggtitle(sprintf("UBC: res=%1.1f",curRes)); 
    outFile <- sprintf("%s/UBC_DimPlot_res%1.1f.%s.pdf",outDir,curRes,dt)
    ggsave(p,file=outFile)

    p <- DotPlot(srat, UBCmarkers, assay="SCT",group.by=groupBy) + 
    ggtitle(ttl)
    outFile <- sprintf("%s/UBC_DotPlot_res%1.1f.%s.pdf",outDir,curRes,dt)
    ggsave(p,file=outFile,width=13,height=6,unit="in")

    cellcat <- unique(mdata[,groupBy])
    plist <- list()
    for (k in cellcat) {
        srat$cur_cluster <- mdata[,groupBy] == k
        table(srat$cur_cluster)
        plist[[k]] <- DimPlot(srat, reduction="umap", 
            cols=c("grey90","#aa0000"),
            group.by="cur_cluster",
            order=TRUE,
            label=FALSE) + ggtitle(k) + 
            xlab("UMAP 1") + ylab("UMAP 2") + 
            theme(plot.title=element_text(size=9)) + NoLegend()
    }
    p <- ggarrange(plotlist=plist)
    ggsave(p,
        file=sprintf("%s/UBC_DimPlot_Single_res%1.1f.%s.pdf",
            outDir,curRes,dt),
        width=8,height=8,unit="in"
    )

}

p <- DimPlot(srat, group.by="sex")
p <- p + ggtitle(sprintf("%s: sex",ttl))
ggsave(p,file=sprintf("%s/UBC_DimPlot_sex_%s.png",outDir,dt))

p <- DimPlot(srat, group.by="dataset_name")
p <- p + ggtitle(sprintf("%s: dataset name",ttl))
ggsave(p,file=sprintf("%s/UBC_DimPlot_datasetname_%s.png",outDir,dt))

outFile <- sprintf("%s/UBC_withClusterAssignments_%s.qs",outDir,dt)
qs::qsave(srat, file=outFile)
