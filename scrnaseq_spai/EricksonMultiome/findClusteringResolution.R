#rm(list=ls())


require(Seurat)
library(clustree)
mb <- "/home/rstudio/isilon/private/projects/MB_multiome/output/QC/20250206/fastmnn/mb_fastmnn.qs"

outDir <- sprintf("/home/rstudio/isilon/private/projects/MB_multiome/output/clustering")
###dt <- format(Sys.Date(),"%y%m%d")
###outDir <- file.path(outDir,dt)

clusterFile <- sprintf("%s/250502/clustering_resolutions.txt", outDir)
outDir <- dirname(clusterFile)

if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE)
}

 cat("Loading Erickson MB data - takes 13 min on Pai lab server\n")
    t0 <- Sys.time()
    #srat <- qs::qread(mb)
    print(Sys.time()-t0)
    cat("Normalizing data\n")
    t0 <- Sys.time()

    full_srat <- srat 

clRes <- c(0.9,1.0)#seq(0.9,1.9,0.2)
  # Perform clustering at multiple resolutions
    for (i in clRes) {
        message("Generating clustering at ", i, " resolution.")
        t0 <- Sys.time()
        srat <- FindClusters(srat, resolution = i)
        print(Sys.time()-t0)
        pattern <- paste0(".*_snn_res.", i)
        matching_column <- grep(pattern, colnames(srat@meta.data), value = TRUE)
        new_name <- paste0("snn_res.", i)
        srat@meta.data[[new_name]] <- srat@meta.data[[matching_column]]
        srat@meta.data[[matching_column]] <- NULL
    }
#srat <- FindClusters(srat, resolution = clRes);

plotClustersSingly <- function(ds, group.by){
plotlist <- list()
md <- ds[[]]
for (cur in unique(md[,group.by])) {
    cat("\tPlotting: ",cur,"\n")
    p <- DimPlot(ds, group.by=group.by, 
        cells.highlight=which(ds[[group.by]]==cur),
        cols.highlight="darkred",cols="#eeeeee", raster=TRUE)
    p <- p + theme(legend.position="none") +
        ggtitle(cur) + 
        theme(plot.title = element_text(hjust = 0.5))
    plotlist[[cur]] <- p
    #ggsave(p,file=sprintf("%s/DimPlot_RLdev_SingleR_%s.pdf",outDir,cur))
}
plotlist
}

cat("Plotting clusters\n")
mdata <- srat[[]]
cols <- colnames(mdata)[grep("^snn_res", colnames(mdata),ignore.case=FALSE)]
cols <- setdiff(cols, paste("snn_res.",c(0.5,0.7,0.8,seq(0.9,1.9,0.2)),sep=""))

for (i in cols) {
    cat("Adding ",i," to metadata\n")
    srat[[i]] <- factor(cl[,i])
}

for (curCol in cols){
   # srat[[curCol]] <- as.character(srat[[curCol]])
    cat(sprintf("Plotting %s\n", curCol))   
    t0  <- Sys.time()
    plotlist <- plotClustersSingly(srat, curCol)
    print(Sys.time()-t0)

    p2 <- ggarrange(plotlist = plotlist)
    ggsave(p2, file=sprintf("%s/DimPlot.%s.png",outDir,curCol), 
        width=30, height=20)
    cat("Done!\n")
}

# p <- DimPlot(srat, group.by="snn_res.1")
# ggsave(p,file=sprintf("%s/DimPlot_snn_res.1.png",outDir))

mdata <- srat[[]]
write.table(mdata,file=sprintf("%s/clustering_assignments.txt",outDir),
    sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)