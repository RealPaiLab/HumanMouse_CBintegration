# explore tumour CNV calls from Akdes
rm(list = ls())


library(ComplexHeatmap)

inFile <- "/.mounts/labs/pailab/private/projects/MB_multiome/input/AkdesHarmanci/finalChrMat_morestringent_v2.rda"
mdataFile <- "/.mounts/labs/pailab/private/projects/MB_multiome/input/AndersErickson/20241212_multiome_metadata.tsv"

outRoot <- "/.mounts/labs/pailab/private/projects/MB_multiome/output"
clusterFile <- sprintf("%s/clustering/250502/clustering_assignments.txt", outRoot)

outDir <- sprintf("%s/CNV",outRoot)
dt <- format(Sys.Date(),"%y%m%d")
outDir <- file.path(outDir,dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE)
}

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


logFile <- sprintf("%s/CNVexplore.log", outDir)
sink(logFile, split=TRUE)
tryCatch({

t0 <- Sys.time()
load(inFile)
print(Sys.time() - t0)

x <- finalChrMat_morestringent
cat("Have %i cells and %i chrom segments\n,",ncol(x),nrow(x))
upos <- regexpr("_",colnames(x))
samps <- substr(colnames(x),1,upos-1)
cat(sprintf("Have %i unique samples\n",length(unique(samps))))

cat("keep only MB samples\n")
mdata <- read.delim(mdataFile, header=TRUE, sep="\t")
mdata <- subset(mdata, methyl_dx %in% c("MB, G4", "MB, G3","MB, SHH CHL AD"))
mdata <- subset(mdata, unique_id %in% samps)

idx <- which(samps %in% mdata$unique_id)
x2 <- x[,idx]
cat(sprintf("%i cells are MB\n", ncol(x2)))

cat("now adding tumour type info\n")
upos <- regexpr("_",colnames(x2))
samps <- substr(colnames(x2),1,upos-1)
midx <- match(samps,mdata$unique_id)
if (all.equal(mdata$unique_id[midx],samps) != TRUE) {
    cat("Row names of metadata and CNV matrix do not match. Fix this.\n")
    browser()
}
ttypes <- mdata[midx,]
ttypes$Cell <- colnames(x2)

cat("Now adding cluster assignments to main Seurat\n")
cl <- read.table(clusterFile, header=TRUE,sep="\t")
y <- as.data.frame(cl[,"snn_res.1"])
y$Cell <- rownames(cl)

x2 <- subset(x2, select = colnames(x2) %in% y$Cell)
y <- subset(y, Cell %in% colnames(x2))

colnames(y)[1] <- "snn_res.1"
midx <- match(colnames(x2), y$Cell)
if (all.equal(colnames(x2), y$Cell[midx]) != TRUE) {
    cat("Column names of CNV matrix and clustering do not match. Fix this.\n")
    browser()
}
clusters <- y[midx,]

cat("making merged dframe\n")
merged <- merge(ttypes, clusters, by = "Cell")

midx <- match(colnames(x2), merged$Cell)
if (all.equal(colnames(x2), merged$Cell[midx]) != TRUE) {
    cat("Column names of CNV matrix and merged data do not match. Fix this.\n")
    browser()
}
x2 <- x2[,midx]


# add coloured bars at the side for each cluster
# create a named vector for the clusters
# created a row annotation for the tumour types
tumour_colors <- c(
    "MB, G3" = as.character(colList["G3"]), 
    "MB, G4" = as.character(colList["G4"]), 
    "MB, SHH CHL AD" = as.character(colList["SHH"])
)

n <- length(unique(mdata))

cat("downsampling")
idx <- which(merged$snn_res.1 %in% c(5,7,9, 28,34,38,20)) #sample(1:ncol(x2), 10000)
x2 <- x2[, idx]
merged <- merged[idx, ]

library(RColorBrewer)
n <- length(unique(merged$snn_res.1))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
clust_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
    rownames(qual_col_pals)))
names(clust_cols) <- unique(merged$snn_res.1)

n <- length(unique(merged$unique_id))
div_col_pals <- brewer.pal.info[brewer.pal.info$category == 'div',]
sample_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
    rownames(qual_col_pals)))
names(sample_cols) <- unique(merged$unique_id)

# create a row annotation for merged$methyl_dx
row_anno <- columnAnnotation(
    Cluster = anno_simple(merged$snn_res.1, col = clust_cols),
    Subgroup=anno_simple(merged$methyl_dx, col=tumour_colors),
    Patient = anno_simple(merged$unique_id, col = sample_cols),
    show_legend = TRUE, height= unit(0.8, "in") 
)

lgd_list <- list(
  Legend(labels = names(clust_cols), title = "cluster", 
    legend_gp = gpar(fill = clust_cols)),
  Legend(labels = names(tumour_colors), title = "subgroup", 
    legend_gp = gpar(fill = tumour_colors))
)

ht_opt$message <- FALSE
ht_opt$ROW_ANNO_PADDING = unit(1,"cm")
ht <- Heatmap(x2, 
               name = "CNV",
               column_order = order(merged$snn_res.1),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               show_column_names = FALSE,
               col = c("blue", "white", "red"),
               row_title = "Chromosomal Segments",
               column_title = "CNVs",
               use_raster = TRUE,
               raster_device = "png",
               top_annotation = row_anno)

pdf(sprintf("%s/CNV_heatmap_clusters.pdf", outDir), width = 16, height = 10)
draw(ht, annotation_legend_list = lgd_list, merge_legend = TRUE)
dev.off()

# Which CNVs are enriched in clusters 5,7,9?
cat("Now looking for enriched CNVs in clusters 5,7,9\n")
clusters_of_interest <- c(5, 7, 9)
cnv_enrichment <- list()
for (cluster in clusters_of_interest) {
    cat(sprintf("Processing cluster %i\n", cluster))
    # Get the cells in the cluster
    cells_in_cluster <- merged$Cell[merged$snn_res.1 == cluster]
    
    # Get the CNV data for these cells
    cnv_data <- x2[, colnames(x2) %in% cells_in_cluster]
    
    # Calculate the mean CNV for each segment
    mean_cnv <- rowMeans(cnv_data, na.rm = TRUE)
    
    # Store the results
    cnv_enrichment[[as.character(cluster)]] <- mean_cnv
}
cnv_enrichment_df <- do.call(cbind, cnv_enrichment)
cnv_enrichment_df <- as.data.frame(cnv_enrichment_df)



}, error = function(e) {
    cat("Error: ", e$message, "\n")
}, finally = {
    sink()
})