rm(list=ls())
library(Seurat)
library(ggplot2)

 mb <- "/home/rstudio/isilon/private/projects/MB_multiome/output/QC/20250206/fastmnn/mb_fastmnn.qs"

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/EYS/Erickson"
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outDir,dt)
if (!file.exists(outDir)) dir.create(outDir)

t0 <- Sys.time()
tumour_srat <- qs::qread(mb)
print(Sys.time()-t0)
DefaultAssay(tumour_srat) <- "SCT"

p <- FeaturePlot(tumour_srat, 
    features=c("EYS","CRX","NRL","BRCA1"),
    order=TRUE,
    alpha=0.8,
    pt.size=1.5)
#p <- p + theme(text=element_text(size=24),
#    axis.text=element_text(size=24)) + xlab("UMAP 1") + ylab("UMAP 2")
#pdf(sprintf("%s/VlnPlot_BRCA_%s.pdf",outDir,dt))
pdf(sprintf("%s/FeaturePlot_EYS_%s.pdf",outDir,dt))
print(p)
dev.off()



