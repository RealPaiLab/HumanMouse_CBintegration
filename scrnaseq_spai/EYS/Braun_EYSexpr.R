rm(list=ls())
library(Seurat)
library(rhdf5)
library(ggplot2)
library(ggpubr)

inDir <- "/home/rstudio/isilon/src/neurodev-genomics/scRNAseq/Braun_2023"
inFile <- sprintf("%s/HumanFetalBrainPool.h5",inDir)

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/EYS/Braun_2023"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outDir,dt)
if (!file.exists(outDir)) dir.create(outDir)

logFile <- sprintf("%s/logfile.txt",outDir)
sink(logFile,split=TRUE)

tryCatch({

cat("* Reading Braun data...\n")
# find the sets using something similar to...

# x<- H5file$new()
# gp <- x[["shoji"]]
# write.table(gp,file="groups.txt")
emb <- h5read(file=inFile,name="shoji/Embedding")
ccl <- h5read(file=inFile,name="shoji/CellClass")
age <- h5read(file=inFile,name="shoji/Age")
genes <- h5read(file=inFile,name="shoji/Gene")
sex <- h5read(file=inFile,name="shoji/Sex")

###cat("Fetching BRCA1 expression\n")
###brca1_idx <- grep("BRCA1",genes)
###brca1 <- h5read(file=inFile,name="shoji/Expression",index=list(brca1_idx,NULL))
###
#valid <- h5read(file=inFile, name="shoji/ValidCells")
# all are valid apparently

cat(sprintf("Read %i cells\n\n", ncol(embed)))


# plot UMAP with subset of 100K randomly sampled cells.
set.seed(123)
#idx <- sample(1:ncol(emb),100000,FALSE)
embed2 <- t(emb) # t(emb[,idx])
cl2 <-  ccl #ccl[idx]
cl2 <- factor(cl2, levels=c("Neural crest",
    "Radial glia","Neuronal IPC",       
    "Glioblast", "Neuroblast","Neuron", 
    "Oligo","Placodes","Vascular","Erythrocyte","Fibroblast","Immune"))

x <- data.frame(x=embed2[,1],y=embed2[,2],Class=cl2)
x$Class <- factor(x$Class)
p <- ggplot(x, aes(x=x,y=y)) + geom_point(aes(colour=Class))
p <- p + xlab("Embedding 1") + ylab("Embedding 2") 
p <- p + ggtitle("Braun 2023: Embeddings")
pdf(sprintf("%s/Embeddings.pdf",outDir))
print(p)
dev.off()

cat("done embeddings\n")
browser()

#### plot BRCA1 expression
###p <- ggplot(x, aes(x=x,y=y)) + geom_point(aes(colour=BRCA1))
###p <- p + xlab("Embedding 1") + ylab("Embedding 2") + ggtitle("Braun 2023: BRCA1 expression")
###pdf(sprintf("%s/BRCA1.pdf",outDir))
###print(p)
###dev.off()
###
###p <- ggplot(x,aes(y=BRCA1,x=Class)) 
###p <- p + geom_violin(aes(fill=Class)) 
###p <- p + xlab("Cluster") + ylab("log2(BRCA1+1)") 
###p <- p + ggtitle(sprintf("Braun et al. 2023: %s cells",
###    prettyNum(nrow(x),big.mark=",")))
###p <- p + theme(legend.position="none",
###    text=element_text(size=20))
###pdf(sprintf("%s/Braun_BRCA1_bygroup.pdf",outDir),
###    width=13,height=6)
###print(p)
###dev.off()

plotGeneByClass <- function(geneName, genes, xprFile,coords,cclass,
    subsample=NULL,oDir) {

    idx <- which(genes == geneName)
    if (length(idx)<1) stop(sprintf("%s: gene not found.\n",geneName))
    cat(sprintf("Reading expression for %s\n",geneName))
    t0 <- Sys.time()
    xpr <- h5read(file=xprFile,
        name="shoji/Expression",index=list(idx,NULL))
    print(Sys.time()-t0)
    
    if (!is.null(subsample)){
        cat(sprintf("Subsampling %i cells\n",subsample))
        idx <- sample(1:ncol(coords),subsample,FALSE)
        embed2 <- coords[,idx]
        cl2 <- cclass[idx]
        xpr2 <- log2(xpr[idx]+1)
    } else {
        embed2 <- coords
        cl2 <- cclass
        xpr2 <- log2(xpr[1,]+1)
    }
    cl2 <- factor(cl2, levels=c(
        "Neural crest","Radial glia","Neuronal IPC",       
        "Glioblast", "Neuroblast","Neuron", 
        "Oligo","Placodes","Vascular","Erythrocyte","Fibroblast","Immune")
    )
    
    x <- data.frame(x=embed2[1,],y=embed2[2,],Class=cl2,GENE=xpr2)

    nPal <- rev(RColorBrewer::brewer.pal(n=5,name="Purples"))
    lvls <- levels(cl2)
    cols <- c(rep(nPal[1],2),nPal[2],nPal[3],nPal[4],nPal[5],
        rep("#aaaaaa",6))
    names(cols) <- lvls

    cat("Making violin plot\n")
    p <- ggplot(x,aes(y=GENE,x=Class)) 
    p <- p + geom_violin(aes(fill=Class)) 
    p <- p + scale_fill_manual(values=cols)
    p <- p + xlab("Cluster") + ylab("log2(expression)") 
    p <- p + ggtitle(sprintf("Braun et al. 2023: %s: %s cells",
        geneName, prettyNum(nrow(x),big.mark=",")))
    p <- p + theme(legend.position="none",
        text=element_text(size=20))
    pdf(sprintf("%s/Braun_%s_bygroup.pdf",outDir,geneName),
        width=15,height=6)
    print(p)
    dev.off()

    cat("Making boxplot plot\n")
    p <- ggplot(x,aes(y=GENE,x=Class)) 
    p <- p + geom_boxplot(aes(fill=Class)) 
    p <- p + scale_fill_manual(values=cols)
    p <- p + xlab("Cluster") + ylab("log2(expression)") 
    p <- p + ggtitle(sprintf("Braun et al. 2023: %s: %s cells",
        geneName, prettyNum(nrow(x),big.mark=",")))
    p <- p + theme(legend.position="none",
        text=element_text(size=40),
        axis.text.x=element_text(angle=45,hjust=1),
        )
    pdf(sprintf("%s/Braun_%s_boxplot_bygroup.pdf",outDir,geneName),
        width=15,height=6)
    print(p)
    dev.off()

   
    x <- subset(x, Class %in% c("Neural crest", "Radial glia", "Neuronal IPC",
        "Glioblast","Neuroblast","Neuron"))
    

    cat("Neuronal lineage only\n")
    p <- ggplot(x,aes(y=GENE,x=Class)) 
    p <- p + geom_boxplot(aes(fill=Class),width=0.4,outlier.size=0.6) 
    p <- p + scale_fill_manual(values=cols)
    p <- p + xlab("Cluster") + ylab("log2(BRCA1)") + coord_flip()
    p <- p +  stat_compare_means(
        aes(label=stat(test)),
        method="wilcox.test",
        method.args=list(alternative="two.sided"),
        label.x.npc = "center",
        tip.length = 0.02,#,
        symnum.args = list(
            cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
            symbols = c("***", "***", "**", "*", "ns")),
        comparisons=list(
            c("Neural crest","Neuron"),
            c("Radial glia", "Neuron"),
            c("Neuronal IPC","Neuron"),
            c("Glioblast","Neuron")
        ),
        size=12
    )



    #p <- p + EnvStats::stat_n_text(size=9)
    p <- p + ggtitle(sprintf("Braun et al. 2023: %s: Neuronal lineage: %s cells",
        geneName, prettyNum(nrow(x),big.mark=",")))
    p <- p + theme_minimal()+ theme(legend.position="none",
        text=element_text(size=40),
        axis.text=element_text(size=40),
        axis.text.x=element_text(angle=45,hjust=1)
    )

    fName <- sprintf("%s/Braun_%s_boxplot_Neuronal_bygroup.pdf",outDir,geneName)
    ggsave(p, file=fName, width=18,height=12,units="in")

    cat("Num cells per cluster\n")
    tbl <- as.integer(table(x$Class))
    tbl <- tbl[tbl>0]
    cat(sprintf("Num cells per group [%i - %i], mean=%1.2f",
        min(tbl),max(tbl),mean(tbl)))

    cat("Statistical tests\n")
    totest <- c("Neural crest","Radial glia","Neuronal IPC","Glioblast")
    neur <- which(x$Class == "Neuron")
    for (k in totest){
        idx <- which(x$Class == k)
        tt <- t.test(x$GENE[idx], x$GENE[neur],alternative="greater")
        cat(sprintf("t-test: %s vs %s\n",k,"Neuron"))
        print(tt)
        cat(sprintf("%s: Mean %1.2f , %s: Mean: %1.2f\n",
            k, mean(x$GENE[idx]),"Neuron", mean(x$GENE[neur])))
        cat(sprintf("WMW(%s > Neuron), p < %1.2e\n",  
            k, tt$p.value))
        cat(sprintf("N1= %i cells; Neuron=%i cells\n\n",
            length(idx),length(neur)))
    }

    ##cat("Plotting UMAP\n")
    ##idx <- sample(1:nrow(x),100000,F)
    ##x <- x[idx,]
    ##x$GENEX <- 2^x$GENE
    ##p <- ggplot(x, aes(x=x,y=y)) + geom_point(aes(fill=GENEX),alpha=0.2)
    ##p <- p + scale_color_viridis_b()
    ##p <- p + xlab("Embedding 1") + ylab("Embedding 2")
    ##p <- p + ggtitle(sprintf("Braun 2023: %s",geneName))
##
    ##png(sprintf("%s/Braun_Embedding_%s.png",outDir,geneName))
    ##print(p)
    ##dev.off()
}

###plotGeneByClass("SOX2",genes,xprFile=inFile,coords=emb,cclass=ccl, oDir=outDir)
plotGeneByClass("GFAP",genes,xprFile=inFile,coords=emb,cclass=ccl, oDir=outDir)
#plotGeneByClass("BRCA1",genes,xprFile=inFile,coords=emb,cclass=ccl, oDir=outDir)

plotGeneByClass("EYS",genes,xprFile=inFile,coords=emb,cclass=ccl, oDir=outDir)

},error=function(ex){
    print(ex)
},finally={
    sink(NULL)
})