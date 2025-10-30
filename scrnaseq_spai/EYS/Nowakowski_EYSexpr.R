# plot EYS in Nowakowski dataset
rm(list=ls())
library(reshape2)
library(ggpubr)
library(RColorBrewer)

inDir <- "/home/rstudio/isilon/src/neurodev-genomics/scRNAseq/Nowakowski_2017"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/EYS/Nowakowski_2017"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outDir,dt)

if (!file.exists(outDir)) dir.create(outDir)
logFile <- sprintf("%s/logFile.txt",outDir)
sink(logFile,split=TRUE)
tryCatch({

cat("* Reading in data...\n")
pheno <- read.delim(sprintf("%s/meta.tsv",inDir),sep="\t",h=T,as.is=T)
dat <- read.delim(
        sprintf("%s/exprMatrix.tsv.gz",inDir),sep="\t",h=T,as.is=T)
cat(sprintf("Nowakowski: %i x %i matrix\n",nrow(dat),ncol(dat)))
cat(sprintf("Pheno: %i cells\n\n",nrow(pheno)))

gn <- dat[,1];
dat <- dat[,-1]
fulldat <- dat

pheno <- pheno[-which(pheno$WGCNAcluster == ""),]
pheno$WGCNAcluster[grep("nEN",pheno$WGCNAcluster)] <- "nEN"
pheno$WGCNAcluster[grep("nIN",pheno$WGCNAcluster)] <- "nIN"
pheno$WGCNAcluster[grep("IN",pheno$WGCNAcluster)] <- "IN"
pheno$WGCNAcluster[grep("EN",pheno$WGCNAcluster)] <- "EN"
pheno$WGCNAcluster[grep("MGE",pheno$WGCNAcluster)] <- "MGE"
pheno$WGCNAcluster[which(pheno$WGCNAcluster %in% c("Astrocyte","Endothelial","Microglia","OPC","Mural"))] <- "non-neuronal"

levelOrder <- c(
    "RG-early","RG-div1","RG-div2","oRG","tRG","vRG",
    "IPC-div1","IPC-div2","MGE","EN","IN","Glyc",
    "non-neuronal"
    )

pheno <- subset(pheno, RegionName %in% "Cortex")

cat("FILTER: RegionName= Cortex: %i cells\n\n",nrow(pheno))
cat("Clusters after consolidating:\n")
print(table(pheno$WGCNAcluster))

nPal <- rev(RColorBrewer::brewer.pal(n=5,name="Purples"))

dat <- melt(dat[grep("EYS",gn),])
colnames(dat) <- c("Cell","EYS")

brca1 <- merge(x=dat,y=pheno,by="Cell")

###tokeep <- c(
###        grep("RG",brca1$WGCNAcluster),
###        grep("IPC",brca1$WGCNAcluster),
###        grep("^EN",brca1$WGCNAcluster),
###        grep("^IN",brca1$WGCNAcluster),
###        grep("non-neuronal",)
###)
###brca1 <- brca1[tokeep,]
brca1$WGCNAcluster <- factor(brca1$WGCNAcluster,levels=levelOrder)
if (any(is.na(brca1$WGCNAcluster))) {
    brca1 <- brca1[-which(is.na(brca1$WGCNAcluster)),]
}
cat(sprintf("EYS: %i cells\n",nrow(brca1)))
print(table(brca1$WGCNAcluster,useNA="always"))

brca1$EYS <- log2(brca1$EYS+1)
umap <- read.delim(sprintf("%s/UMAP.coords.tsv.gz",inDir),h=F)
colnames(umap) <- c("Cell","x","y")
brca1 <- merge(x=brca1,y=umap,by="Cell")
mid <- mean(brca1$EYS)

lvls <- levels(brca1$WGCNAcluster)
cols <- c(rep(nPal[1],6),rep(nPal[2],2),rep(nPal[5],3),rep("#aaaaaa",2))
names(cols) <- lvls

cat("Plotting EYS by group - all lineages\n")
p <- ggplot(brca1,aes(y=EYS,x=WGCNAcluster)) 
p <- p + geom_boxplot(aes(fill=WGCNAcluster))
p <- p + xlab("Cluster") + ylab("log2(EYS)") 
p <- p + ggtitle("Nowakowski et al. 2017")
p <- p + theme_minimal() + 
    theme(
        legend.position="none",
        text=element_text(size=20),
        axis.text.x=element_text(angle=45,hjust=1,size=28),
        axis.text.y=element_text(size=28)
    )
p <- p + expand_limits(y=max(brca1$EYS)*1.5)
p <- p + scale_fill_manual(values=cols)
p <- p + EnvStats::stat_n_text(size=9)
pdf(sprintf("%s/Nowakowski_EYS_bygroup.pdf",outDir),
    width=15,height=6)
print(p)
dev.off()

cat("Plotting EYS by group - selected lineages\n")
lvl2 <- c("RG-early","RG-div1","RG-div2","IPC-div1","IPC-div2","EN","IN")
cols2 <- c(rep(nPal[1],3),rep(nPal[2],2),rep(nPal[5],2))
dat <- subset(brca1, WGCNAcluster %in% lvl2)
dat$grouped <- as.character(dat$WGCNAcluster)
dat$grouped[which(dat$WGCNAcluster %in% c("RG-early","RG-div1","RG-div2"))] <- "RG"
dat$grouped[which(dat$WGCNAcluster %in% c("IPC-div1","IPC-div2"))] <- "IPC"
dat$grouped <- factor(dat$grouped,levels=c("RG","IPC","EN","IN"))
p <- ggplot(dat,aes(y=EYS,x=grouped)) 
p <- p + geom_boxplot(aes(fill=grouped)) + coord_flip()
p <- p + xlab("") + ylab("log2(EYS)") 
p <- p + ggtitle("Nowakowski et al. 2017")
p <- p + theme_minimal() + theme(legend.position="none",
    text=element_text(size=50),
    axis.text.x=element_text(angle=45,hjust=1,size=80),
    axis.text.y=element_text(size=60))
p <- p +  stat_compare_means(
    aes(label=stat(test)),
    method="wilcox.test",
    method.args=list(alternative="two.sided"),
    label.x.npc = "center",
    tip.length = 0.02,#,
    symnum.args = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
        symbols = c("***", "***", "**", "*", "ns")),
    comparisons=list(c("RG","EN"),c("RG","IN"),
        c("IPC","EN"),c("IPC","IN")
    ),
    size=20
    
)
#p <- p + EnvStats::stat_n_text(size=14)
p <- p + scale_fill_manual(values=c(nPal[1],nPal[2],rep(nPal[5],2)))
p <- p + scale_y_continuous(breaks=seq(0,10,5))
ggsave(p, 
    file=sprintf("%s/Nowakowski_EYS_bygroup2.pdf",outDir),
    width=13,height=12,unit="in")

browser()

cat("Running T-tests\n")
totest <- c("RG-early","RG-div1","RG-div2","IPC-div1","IPC-div2")
pvals <- c()
for (k in totest){
    curTest <- brca1$EYS[which(brca1$WGCNAcluster %in% k)]
    mature <- brca1$EYS[which(brca1$WGCNAcluster %in% c("EN","IN"))]

    cat(sprintf("%s vs { EN, IN }\n", k))
    cat(sprintf("# cells: %s = %i cells; EN/IN = %i cells\n",
        k, length(curTest), length(mature)))
    cat(sprintf("EYS: Mean(%s) = %1.2f; Mean(EN/IN) = %1.2f \n",
        k, mean(curTest), mean(mature)))
    cur <- t.test(curTest, mature,alternative="two.sided")
    pvals <- c(pvals, cur$p.value)
    cat(sprintf("t-test test: %s > {EN,IN}: p < %1.2f\n",k, cur$p.value))
    print(cur)
    cat("-----------------\n")
    cat("\n")
}

# UMAP
cat("Plotting UMAP\n")
p1 <- ggplot(brca1,aes(x=x,y=y,col=EYS)) + geom_point(alpha=0.7)
p1 <- p1 + scale_color_gradient2(
    midpoint=mid,
    low="#cccccc",mid="#ffeda0",high="#f03b20",
    space="Lab")
p1 <- p1 + xlab("UMAP 1") + ylab("UMAP 2") 
p1 <- p1 + ggtitle("Nowakowski et al. 2017")
p2 <- ggplot(brca1,aes(x=x,y=y,col=WGCNAcluster)) + geom_point()
p2 <- p2 + xlab("UMAP 1") + ylab("UMAP 2")

p <- ggarrange(plotlist=list(p1,p2))
pdf(sprintf("%s/Nowakowski_UMAP.pdf",outDir),
    width=13,height=6)
print(p)
dev.off()

# plot marker genes
markerGenes <- c("HES1","NES","EOMES","NEUROD6","NEUROD1")
dpos <- regexpr("\\|",gn); gn2 <- substr(gn,1,dpos-1)
idx <- which(gn2 %in% markerGenes)
plist <- list()
for (g in markerGenes) {
    cur <- melt(fulldat[which(gn2 == g),])
    colnames(cur) <- c("Cell","GENE")
    x <- merge(x=cur,y=pheno,by="Cell")
    x$WGCNAcluster <- factor(x$WGCNAcluster,levels=levelOrder)
    x$GENE <- log2(x$GENE+1)
    p <- ggplot(x,aes(x=WGCNAcluster,y=GENE)) +
        geom_boxplot(aes(fill=WGCNAcluster))
    p <- p + xlab("Cell cluster") + ylab("log2(Expression)")
    p <- p + ggtitle(sprintf("Developing cortex: Nowakowski: %s",g))
    p <- p + theme(text=element_text(size=16),
        legend.position="none")
    plist[[g]] <- p
}

p <- ggpubr::ggarrange(plotlist=plist,nrow=3,ncol=2)
outFile <- sprintf("%s/Nowakowski_MarkerGenes.pdf",outDir)
ggsave(p,file=outFile,width=18,height=10,unit="in")

},error=function(ex){
    print(ex)
},finally={
    cat("Closing log.\n")
    sink(NULL)
})