rm(list=ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggpubr)

##inFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/20241031_cca_integ.qs"

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/EYS/AldingerSepp_DevHindbrain"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outDir,dt)
if (!file.exists(outDir)) dir.create(outDir)


inFile <- sprintf("/home/rstudio/isilon/private/projects/BRCA1_BrainDev/AldingerSepp_DevHindbrain/AldingerSepp_withCellTypes0.8.qs")

seppMetadata <-  "/home/rstudio/isilon/private/projects/HumanMouseUBC/Sepp2024_metadata.txt"
seppSexInfo <- sprintf("%s/SeppMetadata_Sex.txt",dirname(seppMetadata))

logFile <- sprintf("%s/logfile.txt",outDir)
sink(logFile,split=TRUE)

tryCatch({    
    cat("Reading human hindbrain set")
    t0 <- Sys.time()
    srat <- qs::qread(inFile,nthreads=8L)
    print(Sys.time()-t0)

    mdata <- srat@meta.data
    age <- mdata$age
    sex <- mdata$sex
    sepp_age <- paste(sub(" wpc", "", mdata$Stage),"PCW")
    age[grep("Sepp",mdata$dataset_name)] <- sepp_age[grep("Sepp",mdata$dataset_name)]

    cat("Age range\n")
    print(table(age,useNA="always"))
    cat("Samples\n")
    cat(sprintf("%i cells, %i genes\n", dim(srat)[2], dim(srat)[1]))

    cat("Reading Sepp sex info\n")
    seppSex <- read.delim(seppSexInfo,sep="\t",h=T,as.is=T)
    idx <- grep("Sepp",mdata$dataset_name)
    matched_sex <- seppSex$sex[
        match(mdata$orig.ident[idx],
        seppSex$orig.ident)]
    mdata$sex[idx] <- matched_sex

    cat("\t* Assigning sex to object...\n")
    srat$sex <- mdata$sex

    cat("getting sample id")
    sepp <- read.delim(seppMetadata)
    sepp$cellID <- rownames(sepp)


    ###sepp$TissueID[which(sepp$TissueID %in% "HUM 15293")] <- "HUM 15293 1" # dummy for pattern matching below
    ###spos <- gregexpr(" ",sepp$TissueID)
    ###tmp <- unlist(sapply(1:length(spos), 
    ###    function(x) { substr(sepp$TissueID[x],1,spos[[x]][2]-1)}))

    mdata$cellID <- rownames(mdata)
    sampleID <- mdata[,c("cellID","sample_id","dataset_name")]
    colnames(sampleID) <- c("cellID","TissueID","dataset_name")

    # just aldinger
    sampleID_ald <- subset(sampleID,    dataset_name=="Aldinger_full_cerebellum_human")
    

    sampleID_both <- rbind(sampleID_ald[,1:2], sepp[,c("cellID","TissueID")])

###    srat <- PrepSCTFindMarkers(srat)
###
###    GOOD_CUTOFF <- 95
###    OK_CUTOFF <- 75
###
###    resCols <- colnames(mdata)[grep("snn_res",colnames(mdata))]
###    idx <- grep("celltype",resCols)
###    if (any(idx)) resCols <- resCols[-idx]
###    ctypes <- list()
###
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
###
###    cat("* Assigning cell cluster identity\n")
###    curRes <- "snn_res.0.8"
###    print(curRes)
###    
###    x <- as.data.frame(
###        table(mdata[,c(curRes,"common_cell_name")])
###    )
###    colnames(x)[1] <- "clusterid"
###    x <- x %>% group_by(clusterid) %>% mutate(pct = (Freq / sum(Freq))*100)
###    x <- as.data.frame(x)
###    x$common_cell_name <- as.character(x$common_cell_name)
###
####    browser()
###    ctype_assign <- list()
###  
###    for (k in unique(x[,1])){
###        y <- x[x[,1]==k,]
###        y <- y[order(y$pct,decreasing=TRUE),]
###        y <- y %>% mutate(cumu=cumsum(pct))
###        maxi <- which.min(abs(y$cumu - GOOD_CUTOFF))
###        # which cell type is most represented in this cluster?
###        cur <- cbind(k,y$common_cell_name[1],y$pct[1])
###        cur <- cbind(cur, paste(y$common_cell_name[1:maxi],collapse=","))
###        if (y$pct[1]>=GOOD_CUTOFF) cur <- cbind(cur,"GOOD")
###        else if (y$pct[1]>=OK_CUTOFF) cur <- cbind(cur,"OK")
###        else cur <- cbind(cur, "MIXED")
###
###        ctype_assign[[k]] <- cur
###    }
###
###   tmp <- as.data.frame(do.call("rbind",ctype_assign))
###    colnames(tmp) <- c("cluster","top_cell_type",
###        "top_proportion","cell_types_to95pctile","cluster_purity")
###     
###    ctypes[[curRes]] <- tmp
##
###
###    curtype <- tmp$top_cell_type[match(mdata[,curRes],tmp$cluster)]
###    srat@meta.data[,cname] <- as.character(curtype);
###    srat@meta.data[,cname] <- factor(
###        srat@meta.data[,cname], 
###        levels=lvls, ordered=TRUE
###    )
###

geneSet <- c("MKI67","WLS","EOMES","LMX1A","","RBFOX3","ATOH1","PAX6",
    "RELN","PDGFRA","TNR","TREM2","CLDN5","ITM2A","SKOR2","SOX14")
 p <- DotPlot(srat, geneSet, assay="SCT",group.by="snn_res.0.8_celltype") + ggtitle("Aldinger/Sepp") 
 fName <- sprintf("%s/AldingerSepp_DotPlot.pdf",outDir)
 ggsave(p,file=fName,width=14,height=4,unit="in") 

###    # UMAP of top cluster assignments
    p <- DimPlot(srat, reduction="umap", 
        cols=cell_colours,
        group.by="snn_res.0.8_celltype",
        order=TRUE,raster=TRUE) 
    p <- p + ggtitle(sprintf("Developing hindbrain (Aldinger, Sepp) %s cells",
        prettyNum(dim(srat)[2],big.mark=","))) + 
        xlab("UMAP 1") + ylab("UMAP 2") + 
        theme(plot.title=element_text(size=12), axis.text=element_text(size=12)
        )
    fName <- sprintf("%s/AldingerSepp_DimPlot_CellType.png",outDir)
    ggsave(p,file=fName, width=6,height=4,unit="in")

cat("getting EYS expression values\n")
    xpr <- srat[["RNA"]]$counts
    brca1 <- xpr[which(rownames(xpr)=="EYS"),]
    brca1 <- log2(brca1+1)
    brca1 <- melt(brca1)

    pdf(sprintf("%s/AldingerSepp_EYS_VlnPlot.pdf",outDir))
    print(p)
    dev.off()

    x <- cbind(rownames(srat[[]]), 
        as.character(srat[[]]$snn_res.0.8_celltype),
        as.character(srat[[]]$sex))
    colnames(x) <- c("Cell","Cluster","Sex")
    brca1$Cell <- rownames(brca1)
    colnames(brca1)[1] <- "EYS"
    y <- merge(x=x,y=brca1,by="Cell") 

    y$Cluster <- factor(y$Cluster,levels=c("RL","VZ","UBC",
        "GN","glutamatergic neuron","GABAergic neuron","interneuron",
        "astrocyte","glia", "brainstem","choroid","endothelial",
        "microglia","oligodendrocyte/OPC"))

    nPal <- rev(RColorBrewer::brewer.pal(n=5,name="Purples"))
    lvls <- levels(y$Cluster)
    cols <- c(rep(nPal[1],2),rep(nPal[2],1),rep(nPal[5],4),
        rep("#aaaaaa",7))
    names(cols) <- lvls

cat("plotting EYS: all clusters\n")
    p <- ggplot(y, aes(x=Cluster,y=EYS)) 
    p <- p + ylab("log2(expression)")
    p <- p + geom_boxplot(aes(fill=Cluster))
    p <- p + theme(axis.text=element_text(size=14), 
        legend.position="none")
    p <- p + scale_fill_manual(values=cols)
    p <- p + ggtitle(sprintf("Aldinger/Sepp: EYS: %s cells",
        prettyNum(nrow(y),big.mark=",")))
    fName <- sprintf("%s/AldingerSepp_EYS_boxplot.pdf",outDir)
    ggsave(p,file=fName,width=16,height=2.5,unit="in")

    cat("Plotting EYS: Neuronal lineage only\n")
    y <- subset(y, Cluster %in% c("RL","UBC","GN","glutamatergic neuron",
        "GABAergic neuron"))
    
    p <- ggplot(y, aes(x=Cluster,y=EYS)) 
    p <- p + ylab("log2(EYS)") + xlab("")  + coord_flip()
    p <- p + ggtitle(sprintf("Aldinger/Sepp: EYS: %s cells",
        prettyNum(nrow(y),big.mark=",")))
    p <- p + geom_boxplot(aes(fill=Cluster),width=0.3,outlier.size=0.7)
    p <- p + theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle=45,hjust=1),
        legend.position="none")
    p <- p + scale_fill_manual(values=cols)
    p <- p + scale_x_discrete(labels=c("glutamatergic neuron"="Glutam.\nneuron",
        "GABAergic neuron"="GABA\nneuron"))
    p <- p + theme_minimal() + theme(
        legend.position="none",
        text=element_text(size=40),
        axis.text.x=element_text(angle=45,hjust=1,size=40),
        axis.text.y=element_text(size=40))
        p <- p +  stat_compare_means(
        aes(label=stat(test)),
        method="wilcox.test",
        method.args=list(alternative="greater"),
        label.x.npc = "center",
        tip.length = 0.02,#,
        symnum.args = list(
            cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
            symbols = c("***", "***", "**", "*", "ns")),
        comparisons=list(
            c("RL","UBC"),c("RL","GN"),
            c("RL","glutamatergic neuron"),c("RL","GABAergic neuron")
        ),
    size=12
)
#p <- p + EnvStats::stat_n_text(size=13)
    fName <- sprintf("%s/AldingerSepp_Neuronal_EYS_boxplot.pdf",outDir)
    ggsave(p,file=fName,width=15,height=8,unit="in")

cat("About to plot sex-specific expression\n")
tmp <- merge(x=y,y=sampleID_both,by.x="Cell",by.y="cellID")
tmp <- tmp[!duplicated(tmp[,c("Sex","TissueID")]),]
cat("Num male samples = ", length(unique(tmp$TissueID[tmp$Sex=="M"])),"\n")
cat("Num female samples = ", length(unique(tmp$TissueID[tmp$Sex=="F"])),"\n")

    p2 <- p+ facet_wrap(~Sex)
    ggsave(p2,
        file=sprintf("%s/AldingerSepp_Neuronal_EYS_boxplot_Sex.pdf",outDir),width=15,height=10,unit="in"
    )
    cat("Cell count by sex and cluster")
    print(table(y$Cluster,y$Sex))

cat("Num cells per cluster\n")
tbl <- as.integer(table(y$Cluster))
print(table(y$Cluster))
tbl <- tbl[tbl>0]
cat(sprintf("Num cells per group [%i - %i], mean=%1.2f",
    min(tbl),max(tbl),mean(tbl)))

    cat("Running T-tests\n")
totest <- c("RL")
pvals <- c()
brca1  <- y
for (k in totest){
    curTest <- brca1$EYS[which(brca1$Cluster %in% k)]
    mature <- brca1$EYS[which(brca1$Cluster %in% c("UBC","GN","glutamatergic neuron","GABAergic neuron"))]

    cat(sprintf("%s vs { UBC, GN, GLU, GABA }\n", k))
    cat(sprintf("# cells: %s = %i cells; mature = %i cells\n",
        k, length(curTest), length(mature)))
    cat(sprintf("EYS: Mean(%s) = %1.2f; Mean(EN/IN) = %1.2f \n",
        k, mean(curTest), mean(mature)))
    cur <- t.test(curTest, mature,alternative="greater")
    pvals <- c(pvals, cur$p.value)
    cat(sprintf("t-test test: %s > {EN,IN}: p < %1.2f\n",k, cur$p.value))
    print(cur)
    cat("-----------------\n")
    cat("\n")
}

  

# UMAP
cat("Plotting UMAP\n")
##embeds <- as.data.frame(srat[["umap"]]@cell.embeddings)
##embeds$Cell <- rownames(embeds)
##brca1 <- xpr[which(rownames(xpr)=="EYS"),]
##brca1 <- log2(brca1+1)
##brca1 <- melt(brca1)
##brca1$Cell<- rownames(brca1)
##colnames(brca1)[1] <- "EYS"
##x <- merge(x=embeds,y=brca1,by="Cell")
##mid <- mean(brca1$EYS)

##p1 <- ggplot(x,aes(x=UMAP_1,y=UMAP_2,fill=EYS)) + 
##    geom_hex(bins=200) #geom_point(alpha=0.7,size=0.5)
##p1 <- p1 + scale_fill_gradient2(
##    midpoint=mid,
##    low="#cccccc",mid="#ffeda0",high="#f03b20",
##    space="Lab")
##p1 <- p1 + xlab("UMAP 1") + ylab("UMAP 2") 
##p1 <- p1 + ggtitle("AldingerSepp")
##fName <- sprintf("%s/AldingerSepp_EYS_FeaturePlot.pdf",outDir)
##ggsave(p1,file=fName,width=6,height=4,unit="in")

###p1 <- FeaturePlot(sratneu,"NES",cols=c("#cccccc","#ffeda0","#f03b20"), #viridis(256),
###    max.cutoff='q98')
###fName <- sprintf("%s/AldingerSepp_NES_FeaturePlot.pdf",outDir)
###ggsave(p1,file=fName,width=5,height=4,unit="in")

###p1 <- FeaturePlot(sratneu,"EYS",cols=c("#cccccc","#ffeda0","#f03b20"), #viridis(256),
###    max.cutoff='q98')
###fName <- sprintf("%s/AldingerSepp_EYS_FeaturePlot.pdf",outDir)
###ggsave(p1,file=fName,width=5,height=4,unit="in")
###
###    cat("Neuronal lineage only\n")
###    sratneu <- subset(srat, snn_res.0.8_celltype %in% 
###        c("RL","UBC","VZ","GN","glutamatergic neuron",
###        "GABAergic neuron","Purkinje cell")
###    )
###sratneu <- RunPCA(sratneu, 
###    features = VariableFeatures(sratneu), npcs = 25)
###sratneu <- RunUMAP(sratneu, reduction = "pca", dims = 1:25)
###
### p <- DimPlot(sratneu, reduction="umap", 
###        cols=cell_colours,
###        group.by="snn_res.0.8_celltype",
###        order=TRUE,raster=TRUE) 
###    p <- p + ggtitle(sprintf("Developing hindbrain (Aldinger, Sepp) %s cells",
###        prettyNum(dim(srat)[2],big.mark=","))) + 
###        xlab("UMAP 1") + ylab("UMAP 2") + 
###        theme(plot.title=element_text(size=12), axis.text=element_text(size=12)
###        )
###    fName <- sprintf("%s/AldingerSepp_DimPlot_CellType.pdf",outDir)
###    ggsave(p,file=fName, width=6,height=4,unit="in")


},error=function(ex){
    print(ex)
},finally={
    sink(NULL)
})