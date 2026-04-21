# Get mouse-only RL integration, subset UBCs, and plot clusters.
# Run DEG analysis and plot heatmap of top markers in each cluster.

library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(harmony)
source("clusterUtils.R")

qsFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_MmRL/20260304/cca/20260304_cca_integ.qs"
gencodeFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/anno/gencode.vM38.basic.annotation.gtf"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_MmRL/UBCclusters"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE)
}

logFile <- sprintf("%s/%s_UBCclusters.log", outDir, dt)
sink(logFile,split=TRUE)

tryCatch({
    cat("reading mouse RL file\n")
    t0 <- Sys.time()
    srat <- qs::qread(qsFile)
    print(Sys.time() - t0)
    # print num cells and num genes
    cat(sprintf("number of cells: %i\n", ncol(srat)))
    cat(sprintf("number of genes: %i\n", nrow(srat)))

    markerGenes <- c("WLS","SOX2","MKI67","EOMES","PAX6","RBFOX3","ITPR1","CA8")

    for (cl in seq(0.2,0.4,0.2)) {
        clustName <- sprintf("snn_res.%s", cl)
        cat(sprintf("plotting clusters for %s\n", clustName))
        p <- DotPlot(srat, features = markerGenes, group.by = clustName) + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(filename = sprintf("%s/MouseRL_%s_dotplot.pdf", outDir, cl),
            plot = p, width = 8, height = 6, dpi = 300)
    }

    
    # keep EOMES+ cells, and then subset. Then cluster just these.
    eomes_plus <- subset(srat, snn_res.0.2 %in% c(6))
    cat(sprintf("number of EOMES+ cells: %i\n", ncol(eomes_plus)))
    # ... (additional subset and clustering steps)
    cat("Running SCTransform and clustering on EOMES+ cells\n")
    t0 <- Sys.time()
    eomes_plus <- SCTransform(eomes_plus, 
        vars.to.regress = c("dataset_name"), verbose = FALSE)
    print(Sys.time() - t0)
    eomes_plus <- RunPCA(eomes_plus, verbose = FALSE, dims=50)
    eb <- ElbowPlot(eomes_plus, ndims=50, reduction="pca")
    ggsave(sprintf("%s/MouseRL_EOMESplus_PCA_elbow.pdf", outDir), plot = eb, width = 6, height = 4, dpi = 300)

    # run Harmony for clustering
    cat("Running Harmony for batch correction\n")
    t0 <- Sys.time()
    harmony <- RunHarmony(eomes_plus, group.by.vars = "dataset_name", assay.use="SCT")
    print(Sys.time() - t0)

    harmony <- RunUMAP(harmony, reduction="harmony", dims = 1:15)
    harmony <- FindNeighbors(harmony, reduction="harmony", dims = 1:15)
    
    harmony$snn_res.0.2 <- NULL
    harmony$snn_res.0.4 <- NULL
    harmony$snn_res.0.6 <- NULL
    harmony <- FindClusters(harmony, resolution = c(0.2,0.4,0.6))

    cat("Plot UMAPs for each resolution\n")
    for (cl in seq(0.2,0.6,0.2)){
        clustName <- sprintf("SCT_snn_res.%s", cl)
        cat(sprintf("plotting clusters for %s\n", clustName))
        p <- DimPlot(harmony, reduction = "umap", group.by = clustName, label = TRUE, label.size = 10) + 
            ggtitle(sprintf("Mouse UBC clusters, %i cells", ncol(harmony)))
        ggsave(filename = sprintf("%s/UBCclusters_%s_DimPlot.pdf", outDir, cl),
            plot = p, width = 8, height = 6, dpi = 300)
    }
    cat("\n")

# plot dataset_name breakdown for each cluster at resolution 0.2
    clustName <- "SCT_snn_res.0.2"
    p2 <- harmony[[]] %>%
        group_by(SCT_snn_res.0.2, dataset_name) %>%
        summarise(count = n()) %>%
        ggplot(aes(x = SCT_snn_res.0.2, y = count, fill = dataset_name)) +
        geom_bar(stat = "identity", position = "stack") +
        theme(legend.position = "top",
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20)) +
        xlab("SCT_snn_res.0.2") + ylab("Cell count") +
        ggtitle("Dataset composition of mouse UBC clusters") +
        theme(plot.title = element_text(hjust = 0.5, size = 20))
    ggsave(filename = sprintf("%s/UBCclusters_SCT_snn_res.0.2_barplot.pdf", outDir),
        plot = p2, width = 8, height = 6, dpi = 300)

# Now similarly show orig.ident breakdown
p2 <- harmony[[]] %>%
  group_by(SCT_snn_res.0.2, orig.ident) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = SCT_snn_res.0.2, y = count, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  xlab("SCT_snn_res.0.2") + ylab("Cell count") +
  ggtitle("orig.ident composition of mouse UBC clusters") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
ggsave(filename = sprintf("%s/UBCclusters_SCT_snn_res.0.2_origident_barplot.pdf", outDir),
    plot = p2, width = 8, height = 6, dpi = 300)


# Run PrepSCTFindMarkers and FindAllMarkers for the 0.2 resolution clusters, and plot heatmap of top 5 markers per cluster.
cat("Running cluster workup\n")
clRes <- "0.2"#seq(0.1,1,0.1)
clusterDir <- sprintf("%s/cluster%1.2f", outDir, as.numeric(clRes))
if (!dir.exists(clusterDir)) {
  dir.create(clusterDir, recursive=FALSE)
}
clusterWorkup(harmony, outDir=clusterDir, clRes=clRes)

# now plot heatmap for the cluster.
res <- 0.2
cat(sprintf("Resolution: %1.2f\n\n", res ));     
resStr <- sprintf("SCT_snn_res.%1.1f", res)
harmony$seurat_clusters <- harmony[[resStr]];# so we can use as a common column    

inFile <- sprintf("%s/Res_%1.1f/DEG_allClusters_SCT_snn_res.%1.1f.tsv", clusterDir, res, res)
de_genes <- read.table(inFile, header=TRUE, sep="\t")    
fcThresh <- 1.5
showTop <- 50

# fetch gencode protein coding genes for annotation using BioC
message("reading gencode")
gencode <- rtracklayer::import(gencodeFile)
gencode <- as.data.frame(gencode)
gencode <- subset(gencode,
    type %in% "gene" & gene_type %in% "protein_coding")

hm <- plot_heatmap(harmony, de_genes, 
    fcThresh=fcThresh, 
    subsetGenes = gencode$gene_name,
    QValThresh=0.05, 
    showTop=showTop,
    genesToInclude=c("SOX4","SOX11"))


# save heatmap
outFile <- sprintf("%s/Res_%1.1f/cluster%1.2f_heatmap_fc%1.2f_top%i.pdf", clusterDir, res, res, fcThresh, showTop)
cat("Saving heatmap to ", outFile, "\n")
pdf(outFile, width=14, height=20)
draw(hm, heatmap_legend_side="right")
dev.off()

DefaultAssay(harmony) <- "RNA"
p <- FeaturePlot(harmony, features = c("EOMES","SOX4","SOX11","LMX1A","OTX2","EYS"), order = TRUE, cols = c("lightgray", "red"))
ggsave(filename = sprintf("%s/EOMESplus_features.pdf", outDir), plot = p, width = 10, height = 10, dpi = 300)


}, error=function(ex){
    print(ex)
}, finally={
    sink()
})



