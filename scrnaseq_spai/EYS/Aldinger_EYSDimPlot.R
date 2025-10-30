rm(list=ls())
library(Seurat)

inFile <- "/home/rstudio/isilon/src/medulloblastoma-genomics/scRNAseq/Aldinger2021_fromLiam/glutamatergic_dev_Liam.RDS"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/EYS/Aldinger"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outDir,dt)
if (!file.exists(outDir)) dir.create(outDir)

liam_srat <- readRDS(inFile)

srat <- liam_srat
srat <- PercentageFeatureSet(srat, pattern="^MT-",col.name="percent.mt")
srat <- SCTransform(srat, vars.to.regress="percent.mt",verbose=FALSE)
cat("\n\nRunPCA and RunUMAP ...");
srat <- RunPCA(srat, verbose = FALSE)
srat <- RunUMAP(srat, dims = 1:20, verbose = FALSE)
		
cat("\n\nFindNeighbors and FindClusters ...");
srat <- FindNeighbors(srat, dims = 1:20, verbose = FALSE)
srat <- FindClusters(srat, verbose = FALSE)

# set Liam cell type to be idents
mdata <- srat@meta.data
Idents(srat) <- mdata$new_cell_type

cat("Printing sample information")
numsamp <- length(unique(mdata$sample_id))
cat(sprintf("%i samples\n",numsamp))
cat(sprintf("%i cells\n",length(mdata$sample_id)))
cat("Age distribution:\n")
print(table(mdata$age))


p <- DimPlot(srat,
    cols=c("#2B2B77","#766BB3","#85A1CA","#345F9D","grey80","gray70","gray40"),pt.size=1)
p <- p + theme(text=element_text(size=24),
    axis.text=element_text(size=24)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(p,file=sprintf("%s/DimPlot.pdf",outDir) ,
width=7,height=7,unit="in")
p <- p + theme(legend.position="none")
ggsave(p,file=sprintf("%s/DimPlot_noLegend.pdf",outDir), width=7,height=7,unit="in")

###x <- FindMarkers(
###    srat, 
###    assay="SCT",
###    logfc.threshold = 0, # default 0.25
###    test.use = "wilcox",
###    min.pct = 0.1, # default 0.1
###    ident.1 = "RL-VZ",
###    ident.2 = "RL-SVZ",
###)
###write.table(x,file=sprintf("%s/DEG_RLvz_RLsvz_%s.tsv",outDir,dt),
###    sep="\t",col=TRUE,row=TRUE,quote=F)

p <- FeaturePlot(srat, 
    features="EYS",#c("BRCA1","OTX2","NEUROD1"),
    order=TRUE,
    alpha=0.8,
    pt.size=1.5)
p <- p + ggtitle("EYS expression")
p <- p + theme(text=element_text(size=24),
    axis.text=element_text(size=24)) + xlab("UMAP 1") + ylab("UMAP 2")
#pdf(sprintf("%s/VlnPlot_BRCA_%s.pdf",outDir,dt))
pdf(sprintf("%s/FeaturePlot_EYS_%s.pdf",outDir,dt))
print(p)
dev.off()


p <- VlnPlot(srat, "EYS")
ggsave(p,file=sprintf("%s/ViolinPlot_EYS_%s.pdf",outDir,dt))

p <- DotPlot(srat, features=c("EYS","OTX2","NEUROD1"),
    assay="SCT",
    group.by="new_cell_type")
pdf(sprintf("%s/DotPlot_%s.pdf",outDir,dt))
print(p)
dev.off()
