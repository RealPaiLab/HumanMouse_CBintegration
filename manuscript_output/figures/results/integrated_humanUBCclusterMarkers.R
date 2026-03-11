# get markers for all human UBC clusters from the human-mouse integrated dataset. 
rm(list=ls())
library(ggplot2)
library(dplyr)

de_genes <- "/home/rstudio/isilon/private/llau/results/integrated/20240715/all_tested_genes.csv"
ranked_markers <- readRDS("/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20241029/cluster_marker_ranking.rds")

in_dir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/diffExpr"
humanOnlyDEG <- sprintf("%s/250424/de_genes.tsv", in_dir)

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/IntegratedUBClusterMarkers"
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive=FALSE)
}
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outDir,dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive=FALSE)
}



cat("reading human-only UBC clusters\n")
humanOnlyDEG <- read.table(humanOnlyDEG, header=TRUE, sep="\t")
humanOnlyDEG <- subset(humanOnlyDEG, p_val_adj < 0.05 & avg_log2FC > 0)

markers <- read.csv(de_genes)

topX <- 10
for (nm in names(ranked_markers)) {
    ranked_markers[[nm]]$cluster <- nm
    cat(sprintf("Cluster %s: { %s }\n", nm, ranked_markers[[nm]]$feature[1:topX] %>% paste(collapse=", ")))
}
integrated_ranked <- do.call(rbind, ranked_markers)

# print the top 10 markers for each cluster in integrated_ranked

# compute the pairwise jaccard similarity for each "cluster" in integrated_ranked with each "ubc_subcluster" in humanOnlyDEG
jaccard_similarity <- matrix(0, nrow=length(unique(integrated_ranked$cluster)), ncol=length(unique(humanOnlyDEG$ubc_subcluster)))
rownames(jaccard_similarity) <- unique(integrated_ranked$cluster)
colnames(jaccard_similarity) <- unique(humanOnlyDEG$ubc_subcluster)
for (i in 1:nrow(jaccard_similarity)) {
    for (j in 1:ncol(jaccard_similarity)) {
        cluster_genes <- integrated_ranked %>% filter(cluster == rownames(jaccard_similarity)[i]) %>% pull(feature)
        ubc_genes <- humanOnlyDEG %>% filter(ubc_subcluster == colnames(jaccard_similarity)[j]) %>% pull(gene)
        intersection <- length(intersect(cluster_genes, ubc_genes))
        union <- length(union(cluster_genes, ubc_genes))
        jaccard_similarity[i, j] <- intersection / union
    }
}

# plot the jaccard similarity as a heatmap using the pheatmap library
library(pheatmap)
library(grid)

rownames(jaccard_similarity) <- sub("_","",rownames(jaccard_similarity))
colnames(jaccard_similarity) <- sub("_","",colnames(jaccard_similarity))
p <- pheatmap(jaccard_similarity,
    cluster_rows=FALSE, cluster_cols=FALSE, 
    display_numbers=TRUE, number_format="%.2f", fontsize=24,
    main="Human Only UBC versus integrated UBC clusters"
)
ggsave(p, filename=sprintf("%s/jaccard_similarity_heatmap.pdf", outDir), width=8, height=6)