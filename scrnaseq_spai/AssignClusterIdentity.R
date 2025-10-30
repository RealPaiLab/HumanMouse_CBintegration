# assign cluster identity to integrated dataset
#rm(list=ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

inFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/20241031_cca_integ.qs"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/assignClusterIdentity"

dt <- format(Sys.Date(),"%y%m%d")

#srat <- qs::qread(inFile,nthreads=8L)

# run clustering at some higher resolutions
#srat2 <- FindClusters(srat, resolution = 0.9)
srat$snn_res.0.9 <- srat2$integrated_snn_res.0.9
srat <- PrepSCTFindMarkers(srat)

# look up metadata column with the clustering res
# for each metadata column, get the cluster assignments and 
# run which.max to ascertain which cell forms the majority
# compute a purity index.

mdata <- srat@meta.data
outFile <- sprintf("%s/cca_integ_metadata.txt",outDir)
if (!file.exists(outFile)){
    write.table(mdata,file=outFile,sep="\t",col=TRUE,row=FALSE,quote=FALSE)
}

# Quang Trinh colour palette
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3',
        '#57C3F3', '#476D87',
        '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
        '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
        '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
        '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
        '#968175'
);

#cell_colours <- my36colors[1:length(unique(mdata$common_cell_name))]

# https://stackoverflow.com/questions/27659831/dplyr-apply-function-table-to-each-column-of-a-data-frame

GOOD_CUTOFF <- 95
OK_CUTOFF <- 75

resCols <- colnames(mdata)[grep("snn_res",colnames(mdata))]
idx <- grep("celltype",resCols)
if (any(idx)) resCols <- resCols[-idx]
ctypes <- list()

cell_groups <- list(
        GABA=c("VZ","Purkinje cell","GABAergic neuron"),
        RL_lineage=c("RL","GCP","GN","UBC","glutamatergic neuron"),
        glia=c("astrocyte","glia","oligodendrocyte/OPC","interneuron"),
        ME=c("immune","microglia","endothelial", "pericyte","choroid",
            "red blood cell","meninges"),
        other=c("isthmic nucleus neuron","neuron (other)","brainstem")
)
lvls <- unlist(cell_groups)
cell_colours <- c(
    brewer.pal(3,"Blues"),
    brewer.pal(5,"YlOrRd"),
    brewer.pal(4,"Greens"),
    brewer.pal(7,"Purples"),
    brewer.pal(3,"Greys")
)

geneSet <- c("MKI67","WLS","EOMES","LMX1A","OTX2","RBFOX3","ATOH1","PAX6",
    "RELN","PDGFRA","TNR","TREM2","CLDN5","ITM2A")

for (curRes in resCols){
    print(curRes)
    #cur <- "snn_res_0.8"
    x <- as.data.frame(
        table(mdata[,c(curRes,"common_cell_name")])
    )
    colnames(x)[1] <- "clusterid"
    x <- x %>% group_by(clusterid) %>% mutate(pct = (Freq / sum(Freq))*100)
    x <- as.data.frame(x)
    x$common_cell_name <- as.character(x$common_cell_name)

#    browser()
    ctype_assign <- list()
  
    for (k in unique(x[,1])){
        y <- x[x[,1]==k,]
        y <- y[order(y$pct,decreasing=TRUE),]
        y <- y %>% mutate(cumu=cumsum(pct))
        maxi <- which.min(abs(y$cumu - GOOD_CUTOFF))
        # which cell type is most represented in this cluster?
        cur <- cbind(k,y$common_cell_name[1],y$pct[1])
        cur <- cbind(cur, paste(y$common_cell_name[1:maxi],collapse=","))
        if (y$pct[1]>=GOOD_CUTOFF) cur <- cbind(cur,"GOOD")
        else if (y$pct[1]>=OK_CUTOFF) cur <- cbind(cur,"OK")
        else cur <- cbind(cur, "MIXED")

        ctype_assign[[k]] <- cur
    }

    tmp <- as.data.frame(do.call("rbind",ctype_assign))
    colnames(tmp) <- c("cluster","top_cell_type",
        "top_proportion","cell_types_to95pctile","cluster_purity")
    
    ctypes[[curRes]] <- tmp

    cname <- sprintf("%s_celltype",curRes)

    curtype <- tmp$top_cell_type[match(mdata[,curRes],tmp$cluster)]
    srat@meta.data[,cname] <- as.character(curtype);
    srat@meta.data[,cname] <- factor(
        srat@meta.data[,cname], 
        levels=lvls, ordered=TRUE
    )

    # UMAP of top cluster assignments
    p <- DimPlot(srat, reduction="umap", 
        cols=cell_colours,
        group.by=cname,
        order=TRUE,raster=TRUE) + ggtitle(curRes) + 
        xlab("UMAP 1") + ylab("UMAP 2") + 
        theme(plot.title=element_text(size=9)
    )
    outFile <- sprintf("%s/%s_topCellAssigned_%s.pdf",outDir,curRes,dt)
    ggsave(p,file=outFile)

     # UMAP of top cluster assignments
    p <- DimPlot(srat, reduction="umap", 
        group.by="dataset_name",
        order=TRUE,raster=TRUE) + ggtitle(curRes) + 
        xlab("UMAP 1") + ylab("UMAP 2") + 
        theme(plot.title=element_text(size=9)
    )
    outFile <- sprintf("%s/%s_dataset_%s.pdf",outDir,curRes,dt)
    ggsave(p,file=outFile)


    # store top assignments and the stats
    outFile <- sprintf("%s/%s_topCell_Assigned_stats_%s.txt",outDir,curRes,dt)
    write.table(ctypes[[curRes]],file=outFile,sep="\t",col=T,row=F,quote=F) 

    # show the proportion plot of cell mixtures under the assigned label
    tmp <- srat@meta.data[,c(cname, "common_cell_name")]
    tmp$common_cell_name <- factor(tmp$common_cell_name, levels=lvls)
    colnames(tmp)[1] <- "umbrella"
    p <- ggplot(tmp, aes(fill=common_cell_name,x=umbrella)) + 
        geom_bar(position="fill") + ggtitle(sprintf("%s: Mixture in top cell  assignment", curRes)) + scale_fill_manual(values=cell_colours)
    p <- p + theme(axis.text=element_text(size=6))
    outFile <- sprintf("%s/%s_topCellAssigned_mixtures_%s.png",outDir,curRes,dt)
    ggsave(p,file=outFile,width=15,height=6,unit="in")

    p <- DotPlot(srat, geneSet, assay="SCT",group.by=cname) + ggtitle(curRes)
    outFile <- sprintf("%s/%s_topCellAssigned_DotPlot_%s.pdf",
            outDir,curRes,dt)
    ggsave(p,file=outFile,width=13,height=6,unit="in")

}

logFile <- sprintf("%s/logfile_%s.txt",outDir,dt)
sink(logFile,split=TRUE)
tryCatch({
    fullDat <- dim(srat[["integrated"]])
    cat("Full cerebellum\n")
    cat(sprintf("Integrated assay:\n%i features, %i cells\n", 
        fullDat[1],fullDat[2]))

#    DefaultAssay(srat) <- "RNA"
    fd <- dim(srat[["RNA"]])
    cat(sprintf("RNA assay: %i features, %i cells\n",fd[1],fd[2]))

    cat("Subset for rhombic lip lineage using 0.8 clustering resolution:")
    # subset for just RL lineage cells
    srat_RL <- subset(srat, 
        snn_res.0.8_celltype %in% c("RL","GCP","GN","UBC",
            "endothelial","microglia","oligodendrocyte/OPC")
    )
    fd <- dim(srat_RL[["SCT"]])
    cat(sprintf("SCT assay: %i features, %i cells\n",fd[1],fd[2]))

}, error=function(ex){
    print(ex)
}, finally={
    print("Closing log.\n")
    sink(NULL)
})


outFile <- sprintf("%s/cca_RLlineage_only_%s.qs",outDir,dt)
message("writing to file")
qs::qsave(srat_RL,file=outFile)