# Run SingleR on Erickson MB clusters

#rm(list=ls())

require(Seurat)
require(ggplot2)
require(ggpubr)
require(SingleR)
require(SingleCellExperiment)
require(scran)
library(dplyr)


###mb <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/tumour/MB_multiome/20250206/fastmnn/mb_fastmnn.qs"
mb <- "/home/rstudio/isilon/private/projects/MB_multiome/output/QC/20250206/fastmnn/mb_fastmnn.qs"
outRoot <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/SingleR_UBCclusters_Erickson"

clusterFile <- "/home/rstudio/isilon/private/projects/MB_multiome/output/clustering/250502/clustering_assignments.txt"

#devFile <- "/home/rstudio/isilon/src/medulloblastoma-genomics/scRNAseq/Aldinger2021_fromLiam/glutamatergic_dev_Liam.RDS"

devFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/assignClusterIdentity/cca_RLlineage_only_241107.qs"

UBC_Seurat <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"

# colour scheme
# WNT: #92c5de
# SHH: #0571b0
# Group 3: #f4a582
# Group 4 #ca0020
colList <- c(
    "WNT"="#92c5de",  
    "SHH"="#0571b0",
    "G3"="#f4a582",
    "G4"="#ca0020",
    "G4/G3"="#ee6211"
)

cellOrder <- c("RL",clrs$Cluster, "UBC",
            "GN","endothelial","microglia",
            "oligodendrocyte/OPC")

UBCcolorF <- "../humanUBCclusters/UBCcolours.txt"
clrs <- read.delim(UBCcolorF,header=TRUE,sep="\t")
clrs[,1] <- paste("UBC",clrs[,1],sep="")
# color scheme for SingleR labels
UBC_pal <- clrs$Colour #RColorBrewer::brewer.pal(9,"Set3")

pal <- c("#D0ABE8",UBC_pal, "blue","#d06ca2",
    RColorBrewer::brewer.pal(3,"Greys"))
names(pal) <- cellOrder

singleRFile <- sprintf("%s/250506/Erickson_singleRPred_250506.RData",outRoot) # sprintf("%s/250407/Erickson_singleRPred_250407.RData",outRoot)

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outRoot,dt)

if (!dir.exists(outDir)) {
    dir.create(outDir,recursive=FALSE)
}

get_norm_expr <- function(srat, rerun_norm = TRUE) {
  if (rerun_norm) {
    srat <- NormalizeData(srat, assay = "RNA")
  }
  norm_mat <- GetAssayData(srat, slot = "data", assay = "RNA")
  return(norm_mat)
}

logFile <- sprintf("%s/logFile_%s.txt",outDir,dt)
sink(logFile,split=TRUE)
tryCatch({
    cat("Loading Erickson MB data - takes 13 min on Pai lab server\n")
    t0 <- Sys.time()
  #  tumour_srat <- qs::qread(mb)
    print(Sys.time()-t0)



    print(sprintf("Tumour data: %i samples, %i cells, %i features",
        length(unique(tumour_srat$sample_id)),
        ncol(tumour_srat),
        nrow(tumour_srat)))
    cat("\n")
    
    DefaultAssay(tumour_srat) <- "SCT"

    cat("assigning cluster identities\n")
    cl <- read.table(clusterFile, header=TRUE,sep="\t")
    if (all.equal(rownames(tumour_srat[[]]),rownames(cl))!=TRUE){
        cat("Row names of metadata and clustering do not match. Fix this. \n")
        browser()
    }
    tumour_srat[["snn_res.1"]] <- factor(cl[,"snn_res.1"])
    tumour_srat$seurat_clusters <- factor(cl[,"snn_res.1"])

if (file.exists(singleRFile)){
    cat("SingleR projections found, loading\n")
    load(singleRFile)
} else {
    cat("SingleR projections not found, running SingleR\n")
    cat("Now loading Aldinger/Sepp dev data\n")
    t0 <- Sys.time()
    rl_srat <- qs::qread(devFile)
    print(Sys.time()-t0)
    # count number of cells and genes in rl_srat
    cat("Number of cells and genes in rl_srat\n")
    cat(sprintf("Cells: %i\n", ncol(rl_srat)))
    cat(sprintf("Genes: %i\n", nrow(rl_srat)))
    
    cat("Now loading UBC clusters\n")
    t0 <- Sys.time()
    cat("\n")

    cat("Now loading UBC clusters\n")
    UBC_srat <- readRDS(UBC_Seurat)
    print(Sys.time()-t0)
    cat("Number of cells and genes in UBC_srat\n")
    cat(sprintf("Cells: %i\n", ncol(UBC_srat)))
    cat(sprintf("Genes: %i\n", nrow(UBC_srat)))
    cat("\n")
    
    DefaultAssay(rl_srat) <- "SCT"    
    DefaultAssay(UBC_srat) <- "SCT"
    
    cat("Keeping genes common to RL and tumour dataset\n")
    counts_rl <- GetAssayData(object = rl_srat, slot = "counts")
    counts_tumour <- GetAssayData(object = tumour_srat, slot = "counts")

    # If a count for a gene in a cell is greater than 0, set as TRUE (= 1)
    nonzero_rl <- counts_rl > 0
    nonzero_tumour <- counts_tumour > 0

    # If 1% or more cells are TRUE, keep the gene. Each TRUE value = 1. Taking the sum of all the cells for that gene
    keep_rl_genes <- Matrix::rowSums(nonzero_rl) >= (0.01*length(Cells(rl_srat)))
    keep_tumour_genes <- Matrix::rowSums(nonzero_tumour) >= (0.01*length(Cells(tumour_srat)))

    # Get the genes names that we are keeping
    rl_genes <- rownames(rl_srat)[keep_rl_genes]
    tumour_genes <- rownames(tumour_srat)[keep_tumour_genes]
    intersecting_genes <- intersect(rl_genes, tumour_genes)  

    rl_srat <- subset(rl_srat, features = intersecting_genes)
    cat("Subsetting tumour data to intersecting genes\n")
    t0 <- Sys.time()
    tumour_srat <- subset(tumour_srat, features = intersecting_genes)
    print(Sys.time()-t0)

    cat("Convert tumour data to SCE object\n")
    t0 <- Sys.time()
    tumour_sce <- as.SingleCellExperiment(tumour_srat)
    print(Sys.time()-t0)

    cat("Harmonize cell cluster labels between the UBC set and the RL set")
    cat("UBC labels\n")
    ubc_labels <- as.data.frame(as.character(UBC_srat[[]]$SCT_snn_res.0.5))
    ubc_labels$CellID <- rownames(UBC_srat[[]])
    colnames(ubc_labels) <- c("labels","CellID")
    table(ubc_labels$labels,useNA="always")

    cat("RL labels\n")
    RL_labels <- as.data.frame(as.character(rl_srat[[]]$snn_res.0.8_celltype))
    RL_labels$CellID <- rownames(rl_srat[[]])
    colnames(RL_labels) <- c("labels","CellID")
    table(RL_labels$labels,useNA="always")

    # create a new label column that combines the UBC and RL labels using the CellID
    # merge the two dataframes on CellID
    merged_labels <- merge(RL_labels, ubc_labels, by = "CellID", all.x=TRUE)
    idx <- which(!is.na(merged_labels$labels.y))
    merged_labels$labels.y[idx]<- sprintf("UBC%s",merged_labels$labels.y[idx])
    merged_labels$labels.x[idx] <- merged_labels$labels.y[idx]

    midx <- match(rownames(rl_srat[[]]),merged_labels$CellID)
    if (all.equal(merged_labels$CellID[midx],rownames(rl_srat[[]]))!= TRUE) {
        stop("CellID mismatch")
    }
    rl_srat[[]]$RLwithUBC <- merged_labels$labels.x[midx]

    browser()
    cat("Convert to SCE object\n")
    rl_sce <- as.SingleCellExperiment(rl_srat)

    cat("Running SingleR\n")
    t0 <- Sys.time()
    trained_model <- trainSingleR(
      rl_sce, 
      labels=rl_sce$RLwithUBC, 
      de.method = "wilcox",
      assay.type = "logcounts"
    )
    print(Sys.time()-t0)

    cat("Classifying SingleR\n")
    t0 <- Sys.time()
    singleRPred <- classifySingleR(
      test = tumour_sce,
      trained = trained_model,
      assay.type = "logcounts"
    )
    print(Sys.time()-t0)
    save(singleRPred, file = sprintf("%s/Erickson_singleRPred_%s.RData",outDir,dt))
    cat("SingleR classification complete\n")
}

mdata <- tumour_srat[[]]
mdata2 <- mdata[,c("subtype","seurat_clusters")]
mdata2$CellID <- rownames(mdata2)

singleRPred$CellID <- rownames(singleRPred)
merged <- merge(mdata2, singleRPred, by = "CellID")

cat("Add SingleR labels to tumour_srat\n")
midx <- match(rownames(mdata), merged$CellID)
if (all.equal(merged$CellID[midx],rownames(mdata))!= TRUE) {
    stop("CellID mismatch")
}
tumour_srat$RLdev_SingleR_pruned.labels <- merged$pruned.labels[midx]
tumour_srat$RLdev_SingleR_labels <- merged$labels[midx]

cat("Plot stacked barplot of RLdev_SingleR_label distribution by tumour subtype")

for (xaxis in c("subtype","seurat_clusters")){
    tbl <- tumour_srat[[]][,c(xaxis,"RLdev_SingleR_labels")]
    colnames(tbl)[1] <- "tcat"

    tbl[,2] <- factor(tbl[,2], levels=cellOrder)

    if (xaxis == "subtype") { wd <- 10; ht <- 10; bs <- 24
    } else { wd <- 15; ht <- 8; bs <- 12 }

    p <- ggplot(tbl,
        aes(x=tcat,fill=RLdev_SingleR_labels)) + 
        geom_bar(position="fill") + 
        scale_fill_manual(values=pal) +
        ylab("Proportion") + 
        theme_minimal(base_size=bs) 
    ggsave(p,file=sprintf("%s/SingleR_PredLabels_by%s.pdf",outDir,xaxis),
        width=wd,height=ht)
}

# plot stacked barplot of subtype distribution by seurat_clusters
cat("Plot stacked barplot of subtype distribution by seurat_clusters")
tbl <- tumour_srat[[]][,c("seurat_clusters","subtype")]
tbl$subtype <- factor(tbl$subtype, 
    levels=c("WNT","SHH","G3","G4","G4/G3"))

  p <- ggplot(tbl,
        aes(x=seurat_clusters,fill=subtype)) + 
        geom_bar(position="fill") + 
        scale_fill_manual(values=colList) +
        ylab("Proportion") + 
        theme_minimal(base_size=12) 
    ggsave(p,file=sprintf("%s/SubtypesByClusters.pdf",outDir),
        width=15,height=8)

# order seurat_clusters by proportion of G3 and then G4 cells
cat("Order seurat_clusters by proportion of G3 and then G4 cells")
tbl <- tumour_srat[[]][,c("seurat_clusters","subtype","RLdev_SingleR_labels")]
tbl$subtype <- factor(tbl$subtype, 
    levels=c("WNT","SHH","G3","G4","G4/G3"))
tbl$RLdev_SingleR_labels <- factor(tbl$RLdev_SingleR_labels,
    levels=cellOrder)
tbl[,1] <- as.character(tbl[,1])
# compute proportion of G4 cells per cluster and convert to data.frame
G4desc <- as.data.frame(tbl %>% 
    group_by(seurat_clusters) %>% 
    summarize(G4_prop = sum(subtype=="G4")/n()) %>%
    arrange(desc(G4_prop)))

UBC0desc <- as.data.frame(tbl %>% 
    group_by(seurat_clusters) %>% 
    summarize(UBC0_prop = sum(RLdev_SingleR_labels=="UBC0")/n()))
UBC0desc$seurat_clusters <- factor(UBC0desc$seurat_clusters,
        levels=as.character(G4desc[,1]))

UBC2desc <- as.data.frame(tbl %>% 
    group_by(seurat_clusters) %>% 
    summarize(UBC2_prop = sum(RLdev_SingleR_labels=="UBC2")/n()))
UBC2desc$seurat_clusters <- factor(UBC2desc$seurat_clusters,
        levels=as.character(G4desc[,1]))

tbl$seurat_clusters <- factor(tbl$seurat_clusters,
        levels=as.character(G4desc[,1]))
p <- ggplot(tbl,
        aes(x=seurat_clusters,fill=subtype)) + 
        geom_bar(position="fill") + 
        scale_fill_manual(values=colList) +
        ylab("Proportion") + 
        theme_minimal(base_size=12) 
    ggsave(p,file=sprintf("%s/SubtypesByClusters_ordered.pdf",outDir),
        width=15,height=8)
# now plot the stacked barplot of RLdev_SingleR_labels distribution by seurat_clusters
cat("Plot stacked barplot of RLdev_SingleR_labels distribution by seurat_clusters")
p <- ggplot(tbl,
        aes(x=seurat_clusters,fill=RLdev_SingleR_labels)) + 
        geom_bar(position="fill") + 
        scale_fill_manual(values=pal) +
        ylab("Proportion") + 
        theme_minimal(base_size=12) +
        ggtitle("SingleR labels by clusters - ordered by decreasing %G4")
    ggsave(p,file=sprintf("%s/RLdev_SingleR_ByClusters_ordered.pdf",outDir),
        width=15,height=8)

# plot a barplot of UBC0_prop by seurat_clusters
cat("Plot barplot of UBC0_prop by seurat_clusters")
p <- ggplot(UBC0desc,
        aes(x=seurat_clusters,y=UBC0_prop)) + 
        geom_bar(stat="identity",fill=pal[which(names(pal)=="UBC0")]) +
        ylab("Proportion") + 
        xlab("seurat_clusters") +
        theme_minimal(base_size=24) +
        ggtitle("Proportion of UBC0 cells by clusters")
    ggsave(p,file=sprintf("%s/UBC0_prop_ByClusters.pdf",outDir),
        width=15,height=8)

# plot a barplot of UBC2_prop by seurat_clusters
cat("Plot barplot of UBC2_prop by seurat_clusters")
p <- ggplot(UBC2desc,
        aes(x=seurat_clusters,y=UBC2_prop)) + 
        geom_bar(stat="identity",fill=pal[which(names(pal)=="UBC2")]) +
        ylab("Proportion") + 
        xlab("seurat_clusters") +
        theme_minimal(base_size=24) +
        ggtitle("Proportion of UBC2 cells by clusters")
    ggsave(p,file=sprintf("%s/UBC2_prop_ByClusters.pdf",outDir),
        width=15,height=8)

p <- DimPlot(tumour_srat,label=TRUE, group.by="snn_res.1")
ggsave(p,file=sprintf("%s/DimPlot.pdf",outDir))
###p <- DimPlot(tumour_srat, group.by="subtype")
###ggsave(p,file=sprintf("%s/DimPlot_subtype.pdf",outDir))

# Using tumour_srat, plot a DimPlot for each level in RLdev_SingleR_label such that cells without the label are not plotted and those with the label are dark red. Compile these plots into a list and use ggarrange to arrange them in a grid.
plotlist <- list()
for (cur in cellOrder){
    cat("Plotting RLdev_SingleR_labels: ",cur,"\n")
    p <- DimPlot(tumour_srat, group.by="RLdev_SingleR_labels", 
        cells.highlight=which(tumour_srat$RLdev_SingleR_labels==cur),
        cols.highlight="darkred",cols="#eeeeee")
    p <- p + theme(legend.position="none") +
        ggtitle(cur) + 
        theme(plot.title = element_text(hjust = 0.5))
    plotlist[[cur]] <- p
    #ggsave(p,file=sprintf("%s/DimPlot_RLdev_SingleR_%s.pdf",outDir,cur))
}
cat("plotting composite\n")
p <- ggarrange(plotlist=plotlist)
ggsave(p,file=sprintf("%s/DimPlot_RLdev_SingleR.pdf",outDir), 
    width=15,height=12)

# plot barplot of RLdev_SingleR_labels in decreasing order by proportion of tumour cells that are labelled by it
cat("Plotting barplot of RLdev_SingleR_labels in decreasing order by proportion of tumour cells that are labelled by it\n")
tbl <- as.data.frame(table(tumour_srat$RLdev_SingleR_labels))
tbl$Freq <- tbl$Freq/sum(tbl$Freq) * 100
p <- ggplot(tbl, aes(x=reorder(Var1,-Freq),y=Freq)) + 
    geom_bar(stat="identity",fill="darkred") +
    ylab("% cells") + 
    xlab("RLdev_SingleR_labels") +
    theme_minimal(base_size=24) +
    ggtitle("Proportion of cells labelled by RLdev_SingleR_labels") +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(p,file=sprintf("%s/RLdev_SingleR_labels.pdf",outDir),
        width=15,height=8)  

# create a dataframe showing for each cluster, the proportion of cells labelled by each RLdev_SingleR_label and by each subtype
cat("Creating dataframe showing for each cluster, the proportion of cells labelled by each RLdev_SingleR_label and by each subtype\n")
tbl <- tumour_srat[[]][,c("seurat_clusters","subtype","RLdev_SingleR_labels")]
tbl$subtype <- factor(tbl$subtype, 
    levels=c("WNT","SHH","G3","G4","G4/G3"))
tbl$RLdev_SingleR_labels <- factor(tbl$RLdev_SingleR_labels,
    levels=cellOrder)
tbl[,1] <- as.character(tbl[,1])
# compute proportion of G4 cells per cluster and convert to data.frame
clusterPropCells <- tbl %>% 
    group_by(seurat_clusters) %>% 
    summarize(G4_prop = sum(subtype=="G4")/n(), 
            G3_prop = sum(subtype=="G3")/n(),
            WNT_prop = sum(subtype=="WNT")/n(),
            SHH_prop = sum(subtype=="SHH")/n(),
            G4G3_prop = sum(subtype=="G4/G3")/n(),
            UBC0_prop = sum(RLdev_SingleR_labels=="UBC0")/n(),
            UBC2_prop = sum(RLdev_SingleR_labels=="UBC2")/n(),
            UBC3_prop = sum(RLdev_SingleR_labels=="UBC3")/n(),
            UBC5_prop = sum(RLdev_SingleR_labels=="UBC5")/n(),
            UBC7_prop = sum(RLdev_SingleR_labels=="UBC7")/n(),
            GN_prop = sum(RLdev_SingleR_labels=="GN")/n()
    ) 
clusterPropCells <- as.data.frame(clusterPropCells)

# create a matrix storing the correlation between G4_prop and UBC0_prop, G3_prop and UBC2_prop, G4_prop and UBC3_prop, G4_prop and UBC5_prop, G4_prop and UBC7_prop
corMat <- matrix(NA, nrow=4, ncol=5)
rownames(corMat) <- c("G4_cor_value","G4_cor_pvalue",
    "G3_cor_value","G3_cor_pvalue")
colnames(corMat) <- c("UBC0_prop","UBC2_prop","UBC3_prop","UBC5_prop","UBC7_prop")
for (i in 1:5){
    cur <- colnames(corMat)[i]
    x <- cor.test(clusterPropCells$G4_prop, clusterPropCells[,cur])
    corMat[1,i] <- x$estimate
    corMat[2,i] <- x$p.value

    x <- cor.test(clusterPropCells$G3_prop, clusterPropCells[,cur])
    corMat[3,i] <- x$estimate
    corMat[4,i] <- x$p.value 
}

# now plot a heatmap of the average expression of the UBC genes in the clusters
cat("Computing average expression of all genes\n")

x <- c("MKI67","TOP2A","E2F7","HMGA2","RARB","E2F3",
       "PBX3","PAX2","SLC6A1","ITPR1","SORCS3",
       "SMAD9","PLAGL1","GREM2","LINC00486","CTNT2",
       "EYS","ICQJ-SCHIP1",
       "PEX5L",
       "RFX3","TMEM163", "MGAT4C",
       "OTX2","PRDM6","SNCAIP")

# take intersection of UBC genes and top genes from the Erickson dataset
y <- intersect(x,rownames(srat))

xpr <- AggregateExpression(
  tumour_srat,
  features = y,
  group.by = "snn_res.1",
)[[1]]

xpr <- t(scale(t(xpr)))

# heatmap colours
col_fun <- circlize::colorRamp2(
  breaks = c(-2,0,2),
  colors = c("#FF00FF","#000000", "#FFFF00")
)

# computing singleR proportions
ubcProp <- tumour_srat[[]] %>%
  select(seurat_clusters, RLdev_SingleR_labels) %>%
  table() %>%
  prop.table(margin = 1) #%>%
  #.[ubc_by_age, ]

stype <- tumour_srat[[]] %>%
  select(seurat_clusters, subtype) %>%
  table() %>%
  prop.table(margin = 1) 

cat("Creating top annotation for heatmap\n")
topannot <- HeatmapAnnotation(
   subtype = anno_barplot(
    stype, 
    gp = gpar(fill = colList[colnames(stype)], border="white", lty="solid"),
    bar_width = 1, 
    height = unit(0.33, "in"),
    axis=FALSE
   ),
  UBCsig = anno_barplot(
    ubcProp, 
    gp = gpar(fill = pal[colnames(ubcProp)], border="white", lty="solid"),
    bar_width = 1, 
    height = unit(2.5, "in"),
    axis=FALSE
)
)

lgd_format <- gpar(fontsize=20)
lgd_ttl <- gpar(fontsize=28, fontface = "bold")
cat("creating legend\n")
lgd_list <- list(
  Legend(labels = colnames(ubcProp), title = "RL Cell", 
    legend_gp = gpar(fill = pal[colnames(ubcProp)]), labels_gp = lgd_format, title_gp = lgd_ttl),
  Legend(labels = colnames(stype), title = "MB subtype", 
    legend_gp = gpar(fill = colList[colnames(stype)]), labels_gp = lgd_format, title_gp = lgd_ttl)
)


# plot heatmap
cat("Plotting heatmap\n")
hm <- Heatmap(xpr,
    col=col_fun,
    top_annotation = topannot,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    cluster_rows = TRUE, 
    cluster_columns = TRUE,
    show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 20),
    row_names_gp = gpar(fontsize = 20),
    heatmap_legend_param = list(labels_gp = lgd_format, title_gp = lgd_ttl)
)
pdf(file = file.path(outDir, "heatmap.pdf"),height=24,width=18)
draw(hm, annotation_legend_list = lgd_list, merge_legend = TRUE)
dev.off()





},error=function(e) {
    cat("Error: ",e$message,"\n")
},finally={
    sink()
})