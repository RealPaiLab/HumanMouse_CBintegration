# heatmap of markers of UBC clusters
rm(list=ls())
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)

#DEGresults <- "/home/rstudio/isilon/private/llau/results/integrated/20240715/all_tested_genes.csv"

DEGresults <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/RLlineage_only/241107/UBC_DEG_integrated_snn_res.0.2_241107.txt"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/RLlineage_only/241107/"
setFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/RLlineage_only/241107/UBC_withClusterAssignments_241107.qs"

gencodeFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/gencode.v44.basic.annotation.gtf"

message("reading gencode")
gencode <- rtracklayer::import(gencodeFile)
gencode <- as.data.frame(gencode)
gencode <- subset(gencode,
    type %in% "gene" & gene_type %in% "protein_coding")


clusterColumn <- "integrated_snn_res.0.2"

dt <- format(Sys.Date(),"%y%m%d")

#outDir <- sprintf("%s/%s",outDir,dt)
#if (!file.exists(outDir)) dir.create(outDir)

##message("Loading integrated RL data")
##if (file.exists(normFile)) {
   srat <- qs::qread(setFile) # human only > normalizedata > sctransform run
##} else {
##    # now getting the expression for these
##    srat <- qs::qread(setFile)
##    srat <- subset(srat, subset = species == "human")
##    srat <- NormalizeData(srat)
##    srat <- SCTransform(srat, 
##        vars.to.regress="CC.Difference"
##    )
##    
##    qs::qsave(srat, 
##        file=sprintf("%s/rl_cca_HumanOnly_NormSCTdone_%s.qs",
##        outDir, dt))
##} 

srat <- RunPCA(srat, verbose = FALSE)
srat <- RunUMAP(srat, dims = 1:30, verbose = FALSE)
srat <- FindNeighbors(srat, dims = 1:30, verbose = FALSE)
srat <- FindClusters(srat, verbose = FALSE)
srat_ubc <- ScaleData(srat,assay="SCT",features=rownames(srat))

#ubc_order <- c("UBC_5","UBC_1","UBC_2","UBC_4","UBC_0","UBC_3")
ubc_order <- unique(srat_ubc@meta.data[,clusterColumn])

#srat_ubc$broad_w_ubc_subtypes <- factor(srat_ubc$broad_w_ubc_subtypes,
#    levels=ubc_order)

# create top annotation rows
dset <- srat_ubc@meta.data[,c(clusterColumn,"dataset_name")]
dset[grep("Aldinger",dset[,2]),2] <- "Aldinger"
dset[grep("Sepp",dset[,2]),2] <- "Sepp"
dset <- as.matrix(table(dset))
dset <- t(apply(dset,1,function(x) x/sum(x)))
ages <- as.matrix(table(srat_ubc@meta.data[,c(clusterColumn,"age")]))
ages <- t(apply(ages,1,function(x) x/sum(x)))
sex <- as.matrix(table(srat_ubc@meta.data[,c(clusterColumn,"sex")]))
sex <- t(apply(sex,1,function(x) x/sum(x)))
pal1 <- RColorBrewer::brewer.pal(ncol(ages),"YlOrRd") # age
pal2 <- c("pink","blue") # sex
pal3 <- RColorBrewer::brewer.pal(3,"Dark2") # dataset

# reorder all
ages <- ages[ubc_order,]
sex <- sex[ubc_order,]
dset <- dset[ubc_order,]

topannot <- HeatmapAnnotation(
    age = anno_barplot(ages, 
        gp = gpar(fill = pal1, border="white",lty="solid"),
        bar_width = 1, 
        height = unit(6, "cm"),
        axis=FALSE
    ),
    sex = anno_barplot(sex,
        gp = gpar(fill = pal2,border="white"),
        bar_width = 1,
        height = unit(2,"cm"),
        axis=FALSE
    ),
    dataset = anno_barplot(dset,
        gp = gpar(fill = pal3,border="white"),
        bar_width = 1,
        height = unit(2,"cm"),
        axis=FALSE)
    #annotation_legend_param = dat = list(labels = colnames(dat)))
)
lgd_list <- list(
    Legend(labels=colnames(ages),title="age",legend_gp=gpar(fill=pal1)),
    Legend(labels=colnames(sex),title="sex",legend_gp=gpar(fill=pal2)),
    Legend(labels=colnames(dset),title="dataset",legend_gp=gpar(fill=pal3))
)

orig_dat <- read.delim(DEGresults,sep="\t",h=T,as.is=T)
for (k in unique(orig_dat$cluster)) {
    cur <- subset(orig_dat, cluster == k)
    cur$obs <- -log10(cur$p_val)
    cur$exp <- -log10(runif(nrow(cur)))
    pdf(sprintf("%s/qqplottest.pdf",outDir))
    qqplot(x=cur$exp,y=cur$obs,xlab="-log10(expected)",ylab="-log10(p-value, observed)") 
    abline(0,1,col="red")
    title(sprintf("QQplot UBC cluster %i",k))
    dev.off()

}
orig_dat <- subset(orig_dat, feature %in% gencode$gene_name)
cat(sprintf("Subsetting for protein-coding genes\n"))
for (k in unique(orig_dat$cluster)) {
    cur <- subset(orig_dat, cluster == k)
    cur$obs <- -log10(cur$p_val)
    cur$exp <- -log10(runif(nrow(cur)))
    pdf(sprintf("%s/qqplottest_cluster%i.pdf",outDir,k))
    qqplot(x=cur$exp,y=cur$obs,xlab="-log10(expected)",ylab="-log10(p-value, observed)") 
    abline(0,1,col="red")
    title(sprintf("QQplot UBC cluster %i",k))
    dev.off()

}
# -----------------------------------------------------------------
# Generate heatmaps for different cutoffs

cutoffs <- c(1,1.5,2)
for (cutoff in cutoffs) {
    message("---------------------------------")
    message(sprintf("cutoff = %1.2f", cutoff))

    topGenes <- list()
    topData <- list()
    message("collecting top  genes per UBC cluster")
    for  (cl in ubc_order) {
        dat <- subset(orig_dat,cluster %in% cl)
        dat <- subset(dat, 
            p_val < 0.05
        )
        #dat <- dat[,c(1,2,6,7,14)]
        dat$avgLogFC <- dat$avg_log2FC#(dat[,2]+dat[,4])/2
        dat <- dat[order(dat$avgLogFC,decreasing=TRUE),]
       # message(cl)
        #print(head(dat[,c("avgLogFC","feature")]))
         dat <- subset(dat,
            avgLogFC > cutoff)
        nr <- min(100,nrow(dat))
        if (nr > 0) {
            topGenes[[cl]] <- dat$feature[1:nr]
            topData[[cl]] <- cbind(rep(cutoff,nr),rep(cl,nr),dat[1:nr,])
        }
        
    }
    topData <- do.call("rbind",topData)
    colnames(topData)[1:2] <- c("avgLogFC_cutoff","topGene_cluster")
    write.table(topData,
        file=sprintf("%s/UBC.topGenes.cutoff%1.2f.%s.txt",
            outDir,cutoff,dt),
        sep="\t",col=TRUE,row=FALSE,quote=FALSE
    )

    message("\tplotting heatmap at cellular resolution")
    feat <- c(unlist(topGenes),
            "SOX4","SOX11","FOXP2","EYS","CRX","NRL")
    
    hm <- DoHeatmap(
        srat_ubc,
        features = feat,
        group.by = clusterColumn,
        size=3
    ) + theme(axis.text.y=element_blank())

    hm <- hm + ggtitle(sprintf("UBC markers (Log2FC > %1.2f)",cutoff))
    ggsave(hm, file=sprintf("%s/heatmap_UBCmarkers_cutoff%1.2f.%s.png",
        outDir,cutoff, dt)
    )

    message("plot heatmap at cluster resolution")
    xpr <- as.matrix(AverageExpression(srat_ubc,
        features=feat,
        assay="RNA",
        group.by=clusterColumn
    )[[1]])
    colnames(xpr) <- sub("-","_",colnames(xpr))
    xpr <- xpr[,ubc_order]

    quants <- quantile(xpr,c(0.1,0.95))
    col_fun <- circlize::colorRamp2(quants,c("#000000","#FFFF00"))

    annoGenes <- c("CRX","SOX4","SOX11","EYS","NRL")
    message("Genes to annotate")
    for (nm in names(topGenes)) {
        cur <- topGenes[[nm]]; ln <- length(cur)
        cur <- cur[1:min(length(cur),5)]
        annoGenes <- c(annoGenes, cur)
        message(sprintf("%s: %i in set; selected = { %s }", 
            nm, ln, paste(cur,collapse=",")))
    }

    # add top genes per cluster on the right
    ha <- rowAnnotation(foo = anno_mark(at = match(annoGenes,rownames(xpr)),
        labels = annoGenes))
    hm2 <- Heatmap(xpr, name="AverageExpression",
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            col = col_fun,
            show_row_names = FALSE,
            column_title=sprintf("UBC markers (log2FC > %1.2f), unscaled",
                cutoff), 
        top_annotation = topannot,
        right_annotation = ha
    )

    png(sprintf("%s/heatmap_AveUBCmarkers_cutoff%1.2f.unscaled.%s.png",
        outDir, cutoff,dt),
        height=11,width=7,unit="in",res=72
    )
    draw(hm2, annotation_legend_list = lgd_list)
    dev.off()

    message("now plot scaled")
    unscaled <- xpr
    xpr <- t(scale(t(xpr)))
    col_fun <- circlize::colorRamp2(
        c(-2,0,2),c("#FF00FF","#000000", "#FFFF00")
       #c(0, max(xpr)), c("#000000","#FFFF00")
    )

     # add top genes per cluster on the right
    ha <- rowAnnotation(foo = anno_mark(at = match(annoGenes,rownames(xpr)),
        labels = annoGenes))
    hm3 <- Heatmap(xpr, name="AverageExpression",
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            col = col_fun,
            column_title=sprintf("UBC markers (log2FC > %1.2f), scaled",
                cutoff),
            show_row_names = FALSE,
            top_annotation = topannot,
            right_annotation = ha
    )

    png(sprintf("%s/heatmap_AveUBCmarkers_cutoff%1.2f.scaled.%s.png",
        outDir, cutoff,dt),
        height=11,width=7,unit="in",res=72
    )    
    draw(hm3, annotation_legend_list=lgd_list)
    dev.off()
}