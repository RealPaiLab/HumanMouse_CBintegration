# getting the list of genes in the regulon of the top PySCENIC TFs
rm(list=ls())

library(reshape2)
# PySCENIC result for integrated Aldinger & Sepp human data
rssFile <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/pyscenic/20240704/cell_type_annot/aldinger_sepp_RL.rss.csv"
regulonDir <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/pyscenic/20241017/"

dat <- read.delim(rssFile,nrow=10,sep=",")
topreg <- dat[which(dat$X == "UBC_1"),-1]
topreg <- melt(topreg)
topreg <- topreg[order(topreg[,2],decreasing=TRUE),]

cat(sprintf("Top 5 regulons: { %s }", paste(topreg$variable[1:5],collapse=",")))
topreg[,1] <- sub("\\.","-",topreg[,1])

message("Collect regulon genes")
topreg <- topreg[1:20,]
regulons <- list()
for (k in topreg$variable){
    cur <- read.delim(sprintf("%s/%s.csv",regulonDir,k),sep=",")
    regulons[[k]] <- cur[,1]
}

# now take the genes in each regulon that are driving the signal.
DEGresults <- "/home/rstudio/isilon/private/llau/results/integrated/20240715/all_tested_genes.csv"

dat <- read.delim(DEGresults,sep=",",h=T,as.is=T)
dat <- subset(dat,cluster %in% "UBC_1")
dat <- subset(dat, 
    Aldinger_full_cerebellum_human_p_val < 0.05 & 
    Sepp_full_cerebellum_human_p_val < 0.05
)
message(sprintf("%i genes with Aldinger Q < 0.05 & Sepp Q < 0.05", 
    nrow(dat)))

dat <- dat[,c(1,2,6,7,14)]
dat$avgLogFC <- (dat[,2]+dat[,4])/2
dat$Ald <- dat$Aldinger_full_cerebellum_human_avg_log2FC
dat$Sepp <- dat$Sepp_full_cerebellum_human_avg_log2FC

dat2 <- subset(dat, avgLogFC >0)
message(sprintf("%i of these with avg log2FC > 0", 
    nrow(dat2)))
dat2 <- dat2$feature

for (k in names(regulons)) {
    cur <- intersect(dat2, regulons[[k]])
    cat(sprintf("%s: %i DEG genes = { %s }\n",k,length(cur),
        paste(cur,collapse=",")))
}
