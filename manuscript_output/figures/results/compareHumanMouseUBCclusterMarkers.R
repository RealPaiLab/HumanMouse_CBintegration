rm(list=ls())



mainDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper"
outDir <- sprintf("%s/compareHumanMouseUBCclusterMarkers", mainDir)

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE);  
}
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE);    
}

outDir_root <- outDir

logFile <- sprintf("%s/%s_compareHumanMouseClusters.log", outDir, dt)
sink(logFile,split=TRUE)

tryCatch({

for (mmDName in c("Sepp", "Vladoiu")) {
for (hsDName in c("Sepp", "Aldinger")) {
cat("\n\n***********************************\n")
cat(sprintf("Comparing %s (mouse) with %s (human)\n", mmDName, hsDName))
cat("***********************************\n")
mouseFile<- sprintf("%s/EachMouseDatasetOnly/260307/%s/%s_DEG_results.tsv",mainDir, mmDName, mmDName)
humanFile <- sprintf("%s/EachHumanDatasetOnly/260307/%s/%s_DEG_results.tsv", mainDir, hsDName, hsDName)

outDir <- sprintf("%s/%s_vs_%s", outDir_root, hsDName, mmDName)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE);
}



mm <- read.delim(mouseFile)
hs <- read.delim(humanFile)

# compare markers in pairwise clusters for human and mouse
# compute jaccard

a <-  mm; a_cl <- unique(a$seurat_cluster)
b <- hs; b_cl <- unique(b$seurat_cluster)

a <- subset(a, p_val_adj < 0.1 & avg_log2FC > 0)
b <- subset(b, p_val_adj < 0.1 & avg_log2FC > 0)
cat(sprintf("Genes sig upreg: %i in Mouse (Sepp), %i in Human (Aldinger)\n", nrow(a), nrow(b)))

# compute the overlap in upregulated genes for all pairwise combinations of clusters
common <- list()
mat <- matrix(0, 
    nrow=length(a_cl), ncol=length(b_cl), 
    dimnames=list(paste0("Mouse_cl",a_cl), paste0("Human_cl",b_cl)))
for (i in a_cl) {
    for (j in b_cl) {
        x <- a$gene[which(a$seurat_cluster == i)]
        y <- b$gene[which(b$seurat_cluster == j)]
        
        nm <- sprintf("Mouse_cluster%i_Human_cluster%i", i,j)
        jac <- length(intersect(x,y))/length(union(x,y))
        common[[nm]] <- jac
        mat[which(a_cl == i), which(b_cl == j)] <- jac

        cat(sprintf("Mouse cl %i & Human cl %i: %i and %i genes, %i intersection: %1.2f Jaccard\n",  i, j,
            length(x), length(y), length(intersect(x,y)), jac))

    }
}

# now plot a heatmap of the Jaccard indices for the overlap of upregulated genes between clusters in the two datasets.
# make the text font larger
library(pheatmap)

p <- pheatmap(mat, cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=TRUE, number_format="%.2f", 
    fontsize=18,
    main="Jaccard index of upregulated genes between Aldinger and Sepp clusters")
ggsave(p, file=sprintf("%s/%s_%s_Jaccard_heatmap.pdf", outDir, hsDName, mmDName), width=6, height=5)

# find the pair with the highest jaccard index
top_pair <- names(common)[which.max(unlist(common))]
cat(sprintf("Best pair: %s with Jaccard index %1.2f\n", top_pair, common[[top_pair]]))

# create a table of the genes for each of the top pairs.
df <- list()
for (cur in top_pair) {
    # parse out cluster numbers from top_pairs
    tmp <- sub("Mouse_cluster", "", cur)
    tmp <- sub("Human_cluster", "", tmp)
    upos <- regexpr("_", tmp)
    clA <- as.numeric(substr(tmp, 1, upos-1))
    clB <- as.numeric(substr(tmp, upos+1, nchar(tmp)))
    x <- a$gene[which(a$seurat_cluster == clA)]
    y <- b$gene[which(b$seurat_cluster == clB)]

    common_genes <- intersect(x,y)
    df[[cur]] <- data.frame(gene=common_genes, cluster_pair=cur)

    venn.plot <- draw.pairwise.venn(area1 = length(x), area2 = length(y), cross.area = length(intersect(x,y)), 
        category = c(sprintf("Mouse cl %i", clA), sprintf("Human cl %i", clB)), 
        fill = c("red", "blue"), alpha = 0.5, cat.pos = c(-20, 20), cat.dist = 0.05, scaled = FALSE)

    ggsave(venn.plot, file=sprintf("%s/%s_VennDiagram.pdf", outDir, cur), width=5, height=5)
    write.table(common_genes, file=sprintf("%s/%s_common_genes.txt", outDir, cur), quote=FALSE, row.names=FALSE, col.names=FALSE)
}
}
}


},error=function(ex){
    print(ex)
}, finally={
    sink(NULL)
})