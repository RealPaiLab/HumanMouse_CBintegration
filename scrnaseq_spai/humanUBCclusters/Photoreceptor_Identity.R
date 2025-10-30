# project photoreceptor gene sets onto the UBC clusters and see which clusters express those genes the most.
rm(list=ls())
library(readxl)
library(ggplot2)
library(reshape2)
library(dplyr)

srat_file <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/Photoreceptor"

annoDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/anno"
garancher <- sprintf("%s/Garancher_Photoreceptor_Genes.txt", annoDir)
descartes <- sprintf("%s/suppl_table_5_smith_2022_nature.xlsx", annoDir)

colourFile <- "UBCcolours_250526.txt"
clrs <- read.delim(colourFile,header=TRUE,sep="\t")

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outDir,dt)

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

cat("Reading photoreceptor gene sets\n")
gar <- read.delim(garancher, header=FALSE, sep="\t", as.is=TRUE)
desc <- as.data.frame(read_excel(descartes, sheet="Photoreceptor_gene_set",
    skip=2, col_names=TRUE))

cat("Reading human UBC cluster file\n")
t0 <- Sys.time()
srat <- readRDS(srat_file)
print(Sys.time() - t0)

# remove cells from cluster 6 and cluster 8
cat("Removing clusters 6 and 8 from the UBC object\n")
# remove clusters 6 and 8 from the srat object
srat <- subset(srat, subset = SCT_snn_res.0.5 != 6)
srat <- subset(srat, subset = SCT_snn_res.0.5 != 8)

# set the cluster factor levels to the colours
clrs <- clrs[-which(clrs$Cluster %in% c(6,8)),]
# set the cluster factor levels to the colours
cat("Setting cluster factor levels\n")

# get names of genes in srat
genes <- rownames(srat[["RNA"]]@counts)
gar1 <- intersect(gar[,1], genes)
desc1 <- intersect(as.character(desc[,1]), genes)

cat("Adding photoreceptor scores to the UBC object\n")
srat <- AddModuleScore(srat, features = list(gar1), name = "Garancher", nbin=10)
srat <- AddModuleScore(srat, features = list(desc1), 
    name = "Descartes", nbin=10)

# Plot a violin plot of the Garancher and Descartes scores by cluster
# use the ggplot2 package to create the plots

# Function to plot the scores by cluster
# and save the plots to the output directory
# Note: the scores are stored in the metadata of the srat object
plotPhotoreceptorScores <- function(srat, scoreName) {
  scores <- srat[[scoreName]]
  scores <- as.data.frame(scores)
    scores$Cell <- rownames(scores)

  cl <- srat[["SCT_snn_res.0.5"]]
  cl$Cell <- rownames(cl)

  # Merge scores with clusters
  mergedData <- merge(scores, cl, by = "Cell")
  mergedData$SCT_snn_res.0.5 <- factor(mergedData$SCT_snn_res.0.5, 
                            levels = clrs$Cluster)
  
  # Create the plot
  p <- ggplot(mergedData, aes_string(x = "SCT_snn_res.0.5", y = scoreName)) +
    geom_violin(aes(fill = SCT_snn_res.0.5), alpha = 0.7) +
    scale_fill_manual(values = clrs$Colour) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    labs(title = paste("Photoreceptor Scores by Cluster:", scoreName),
         x = "Cluster",
         y = paste(scoreName, "Score")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
# add pairwise WMW pvalues
    p <- p + stat_compare_means(aes(group = SCT_snn_res.0.5), 
                                label = "p.signif", 
                                method = "wilcox.test", 
                                hide.ns = TRUE, 
                                size = 3.5) 
  return(p)
}

x <- paste("UBC", srat[["SCT_snn_res.0.5"]][,1], sep = "")
clrs$Cluster2 <- paste("UBC", clrs$Cluster, sep = "")
srat$clusters <- factor(x, levels = clrs$Cluster2)

p1 <- plotPhotoreceptorScores(srat, "Garancher1")
# draw a horizontal red line at 0.1
p1 <- p1 + geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
    ggtitle("Garancher Photoreceptor Scores by Cluster")
ggsave(filename = sprintf("%s/Garancher_scores_by_cluster.pdf", outDir), 
       plot = p1, width = 10, height = 6)

p2 <- plotPhotoreceptorScores(srat, "Descartes1")
# draw a horizontal red line at 0.05
p2 <- p2 + geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    ggtitle("Descartes Photoreceptor Scores by Cluster")
ggsave(filename = sprintf("%s/Descartes_scores_by_cluster.pdf", outDir),
       plot = p2, width = 10, height = 6)

# Plot a dotplot of Garancher genes by cluster
dplot <- c("EOMES","SOX2","MKI67","WLS","TBR1", "ITPR1","SMAD9","EYS",gar1)
dplot <- setdiff(dplot, "MAFK")
p <- DotPlot(srat, features = dplot,
    group.by = "clusters", 
    dot.scale = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Garancher Photoreceptor Genes by Cluster")
ggsave(filename = sprintf("%s/Dotplot_Garancher_genes_by_cluster.pdf", outDir), 
       plot = p, width = 10, height = 6)

# plot the violin plot of ages by cluster
age <- as.integer(trimws(sub("PCW","",srat$age)))
srat$age_num <- age

tmp <- srat[["age_num"]]
tmp$Cell <- rownames(tmp)

tmp2 <- srat[["clusters"]]
tmp2$Cell <- rownames(tmp2)
mergedData <- merge(tmp, tmp2, by = "Cell")
mergedData$clusters <- factor(mergedData$clusters, 
                            levels = clrs$Cluster2)
p_age <- ggplot(mergedData, aes(x = clusters, y = age_num)) +
    geom_violin(aes(fill = clusters), alpha = 0.7) +
    scale_fill_manual(values = clrs$Colour) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    labs(title = "Age by Cluster",
         x = "Cluster",
         y = "Age (PCW)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(filename = sprintf("%s/Age_by_cluster.pdf", outDir), 
       plot = p_age, width = 10, height = 6)
