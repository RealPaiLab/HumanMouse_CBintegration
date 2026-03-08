# plot heatmap of UBC markers across subclusters
library(tidyverse)
library(patchwork)
library(Seurat)
library(ComplexHeatmap)

degFile <- "/home/rstudio/isilon/private/llau/results/integrated/20240715/sig_ubc_markers.csv"
featRank <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20241029/cluster_marker_ranking.rds"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/UBCmarkers"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

logFile <- sprintf("%s/heatmap_log.txt", outDir)
sink(logFile, split=TRUE)

tryCatch({

cat("Loading UBC data\n")
srat <- qs::qread("/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20241103/human_ubc_srat_merged.qs")

ald <- subset(ubc_cells, subset = dataset_name == "Aldinger_full_cerebellum_human")
sepp <- subset(ubc_cells, subset = dataset_name == "Sepp_full_cerebellum_human")

cat("reading deg file\n")
deg <- read.delim(degFile, sep = ",", header = TRUE)
feature_ranking <- readRDS(featRank)

ubc_clusters <- levels(ubc_cells$new_clusters)
unique_markers <- map(
  .x = ubc_subclusters,
  .f = \(clust) {
    x <- ranked_markers[[clust]]
    nr <- min(nrow(x), 100)
    x <- pull(x, feature) %>%
      head(nr)
    return(x)
  }
) %>%
  unlist() %>%
  unique()

srat <- ald
cat("Getting average expression for cluster markers\n")
# get data slot, then manually scale the results
xpr <- AverageExpression(
  srat,
  assays = "SCT",
  features = unique_markers,
  group.by = "new_clusters",
  slot = "data"
)[[1]]
colnames(xpr) <- sub("-","_", colnames(xpr))
xpr <- xpr[, ubc_clusters]
xpr <- t(scale(t(xpr)))

# heatmap colours
col_fun <- circlize::colorRamp2(
  breaks = c(-2,0,2),
  colors = c("#FF00FF","#000000", "#FFFF00")
)

#### label top genes
###anno_genes <- map(
###  .x = feature_ranking,
###  .f = \(df) {slice_min(df, order_by = combined_rank, n = 5)}
###)
#### select the feature column from each entry in anno_genes and unlist into a single vector
###anno_genes <- map(anno_genes, "feature") %>%
###  unlist() %>%
###  unname() %>%
###  c("OTX2", "CBFA2T2", "SOX4", "SOX11", "LMX1A", "EOMES")
ha <- rowAnnotation(
  foo = anno_mark(at = 1:nrow(xpr), labels = rownames(xpr),
    labels_gp=gpar(fontsize = 24, col = "black"))
)

# order by cluster
hm <- Heatmap(
  xpr,
  name="scaled\nexpression",
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  col = col_fun,
  show_row_names = FALSE,
  right_annotation = ha,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 20)
)

pdf(
  file = file.path(outDir, "avg_expr_heatmap.pdf"),
  width = 14,
  height = 20
)
draw(hm) # annotation_legend_list = lgd_list, merge_legend = TRUE)
dev.off()

}, error=function(ex){
    print(ex)
}, finally={
    sink(NULL)
})
