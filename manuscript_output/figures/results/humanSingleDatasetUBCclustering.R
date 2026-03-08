
rm(list=ls())
library(Seurat)
library(ggplot2)

aldFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/EachHumanDatasetOnly"
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE);    
} 

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE);
}

logFile <- sprintf("%s/%s_eachhuman_clusters.log", outDir, dt)
sink(logFile,split=TRUE)

tryCatch({
srat <- readRDS(aldFile)
srat_mega <- srat
cat(sprintf("Read UBC data with %d cells and %d genes\n", ncol(srat), nrow(srat)))

# subset for cells with PAX6 and EOMES

degMarkers <- list()
for (dset in c("Aldinger","Sepp")) { #},"Aldinger")) {
    cat(sprintf("\n\n**************** %s\n", dset))
    curDir <- sprintf("%s/%s",outDir,dset)
    if (!file.exists(curDir)) dir.create(curDir,recursive=FALSE)

    cat("Subsetting for dataset, and PAX6+EOMES+ cells\n")
    srat <- subset(srat_mega, dataset_name == sprintf("%s_full_cerebellum_human",dset))
    srat <- subset(srat, subset = PAX6 > 0 & EOMES > 0)
    cat(sprintf("%i cells in %s that are PAX6+ EOMES+\n", ncol(srat), dset))

    cat("Running SCTransform\n")
    srat$percent.mt <- PercentageFeatureSet(srat, pattern = "^MT-")
    srat <- SCTransform(srat, assay="RNA", verbose=FALSE,
        vars.to.regress = c("percent.mt","CC.Difference"),
        return.only.var.genes = FALSE
    )
    cat("Running PCA")
    srat <- RunPCA(srat, assay="SCT", verbose=FALSE)
    eb <- ElbowPlot(srat, ndims=50)
    ggsave(eb, file=sprintf("%s/%s_ElbowPlot.pdf", curDir, dset))

    if (dset == "Sepp") ndims <- 10
    else ndims <- 20
    srat <- FindNeighbors(srat, reduction="pca", dims=1:ndims, verbose=FALSE)
    srat <- RunUMAP(srat, reduction="pca", dims=1:ndims, verbose=FALSE)

    clRes <- 1 #c(0.5, 1, 1.5, 2)
    srat <- FindClusters(srat, resolution=clRes, verbose=FALSE) 
    for (rr in clRes){
        srat$seurat_clusters <- srat[[paste0("SCT_snn_res.",rr)]];# so we can use as a common column
        p <- DimPlot(srat, group.by="seurat_clusters") + xlab("UMAP 1") + ylab("UMAP 2") 
        p <- p + ggtitle(sprintf("%s: Cluster resolution %s", dset, rr))
        p <- p + theme_minimal(base_size = 18)
        ggsave(p, file=sprintf("%s/%s_DimPlot_res%s.pdf", curDir, dset, rr), width=5, height=5)
    }

    clSel <- "SCT_snn_res.1"
    cl2test <- 0
    if (dset == "Sepp"){
        srat$seurat_clusters <- srat[[clSel]] # gives three relatively tight clusters
        cl2test <- 0:2
        cat(sprintf("For Sepp, using clSel = %s and cl2test = {%s}", 
            clSel, paste(cl2test,collapse=",") ))
    } else {
        cat("decide for Aldinger")
        srat$seurat_clusters <- srat[[clSel]] # gives three very loose clusters
        cl2test <- 0:2
    }

    cat("Running DEG analysis\n")
    srat <- PrepSCTFindMarkers(srat, assay = "SCT")
    source("clusterUtils.R")
    Idents(srat) <- clSel # 3 clusters 
    deg <- run_all_de(
      srat = srat,
      clusters = cl2test,
      grouping.var = NULL
    )
    write.table(deg, file=sprintf("%s/%s_DEG_results.tsv", curDir, dset), sep="\t", quote=FALSE, row.names=FALSE)

    cat("Plotting volcano")
    pList <- plotVolcano_usingDEG(deg, showTopGenes=50, logFCcutoff=0.25,titlePfx=dset,base_size=20)
    for (nm in names(pList)) {
        p <- pList[[nm]]
        ggsave(p, file=sprintf("%s/%s_Volcano_%s.pdf",curDir,dset, nm))
    }
    degMarkers[[dset]] <- deg
}

# compute the overlap of significantly upregulated markers in the two datasets. use a log2 fold change cutoff of 0 
# and an adjusted p-value cutoff of 0.05 to define significant upregulation. 
# report the number of markers for each pairwise 
a <- degMarkers[["Aldinger"]]; a_cl <- unique(a$seurat_cluster)
b <- degMarkers[["Sepp"]]; b_cl <- unique(b$seurat_cluster)

a <- subset(a, p_val_adj < 0.1 & avg_log2FC > 0)
b <- subset(b, p_val_adj < 0.1 & avg_log2FC > 0)
cat(sprintf("Genes sig upreg: %i in Aldinger, %i in Sepp\n", nrow(a), nrow(b)))

# compute the overlap in upregulated genes for all pairwise combinations of clusters
common <- list()
mat <- matrix(0, 
    nrow=length(a_cl), ncol=length(b_cl), 
    dimnames=list(paste0("Ald_cl",a_cl), paste0("Sepp_cl",b_cl)))
for (i in a_cl) {
    for (j in b_cl) {
        x <- a$gene[which(a$seurat_cluster == i)]
        y <- b$gene[which(b$seurat_cluster == j)]
        cat(sprintf("Aldinger cl %i & Sepp cl %i: %i and %i genes\n",  i, j,
            length(x), length(y)))

        nm <- sprintf("Ald_cluster%i_Sepp_cluster%i", i,j)
        jac <- length(intersect(x,y))/length(union(x,y))
        common[[nm]] <- jac
        mat[which(a_cl == i), which(b_cl == j)] <- jac
    }
}

# now plot a heatmap of the Jaccard indices for the overlap of upregulated genes between clusters in the two datasets.
# make the text font larger
library(pheatmap)
p <- pheatmap(mat, cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=TRUE, number_format="%.2f", 
    fontsize=18,
    main="Jaccard index of upregulated genes between Aldinger and Sepp clusters")
ggsave(p, file=sprintf("%s/Aldinger_Sepp_Jaccard_heatmap.pdf", outDir), width=6, height=5)

# plot a Venn diagram of the overlap of upregulated genes between the two datasets, using the cluster pair with the highest Jaccard index.
library(VennDiagram)
# use the mat object to get the best pair of clusters with the highest Jaccard index

top_pairs <- c("Ald_cluster0_Sepp_cluster2","Ald_cluster1_Sepp_cluster1","Ald_cluster2_Sepp_cluster0")
# create a table of the genes for each of the top pairs.
df <- list()
for (cur in top_pairs) {
    # parse out cluster numbers from top_pairs
    tmp <- sub("Ald_cluster", "", cur)
    tmp <- sub("Sepp_cluster", "", tmp)
    upos <- regexpr("_", tmp)
    clA <- as.numeric(substr(tmp, 1, upos-1))
    clB <- as.numeric(substr(tmp, upos+1, nchar(tmp)))
    x <- a$gene[which(a$seurat_cluster == clA)]
    y <- b$gene[which(b$seurat_cluster == clB)]

    common_genes <- intersect(x,y)
    df[[cur]] <- data.frame(gene=common_genes, cluster_pair=cur)

    venn.plot <- draw.pairwise.venn(area1 = length(x), area2 = length(y), cross.area = length(intersect(x,y)), 
        category = c(sprintf("Aldinger cl %i", clA), sprintf("Sepp cl %i", clB)), 
        fill = c("red", "blue"), alpha = 0.5, cat.pos = c(-20, 20), cat.dist = 0.05, scaled = FALSE)

    ggsave(venn.plot, file=sprintf("%s/%s_VennDiagram.pdf", outDir, cur), width=5, height=5)
}
df <- do.call(rbind, df)
write.table(df, file=sprintf("%s/common_upregulated_genes_top_ClusterPairs.tsv", outDir), 
    sep="\t", quote=FALSE, row.names=FALSE)

}, error=function(ex){
    print(ex)
}, finally={
    sink()
})