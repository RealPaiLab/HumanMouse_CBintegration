## -----------------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)
library(Seurat)


## -----------------------------------------------------------------------------
thesis_dir <- "../../"
source(file.path(thesis_dir, "./utils.R"))
source(file.path(thesis_dir, "../software/utilities/plotting.R"))

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/UBCmarkers"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}


my_pals <- get_custom_pals()
topX <- 70


## -----------------------------------------------------------------------------
# load Seurat object for human UBCs
human_ubcs <- qs::qread("/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20241103/human_ubc_srat_merged.qs")


## -----------------------------------------------------------------------------
ubc_subclusters <- paste0("UBC_", c(1,5,0,4,3,2))

# load top markers in each cluster across both datasets (Aldinger/Sepp)
# based on Leo's differential expression results from 2024-07-15
ranked_markers <- readRDS("/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20241029/cluster_marker_ranking.rds")
unique_markers <- map(
  .x = ubc_subclusters,
  .f = \(clust) {
    x <- ranked_markers[[clust]]
    nr <- min(nrow(x), topX)
    x <- pull(x, feature) %>%
      head(nr)
    return(x)
  }
) %>%
  unlist() %>%
  unique()

## -----------------------------------------------------------------------------
# heatmap annotations

# create top annotation rows
ages <- human_ubcs[[]] %>%
  select(subclusters, human_age) %>%
  table() %>%
  prop.table(margin = 1) %>%
  .[ubc_subclusters, ]
sex <- human_ubcs[[]] %>%
  select(subclusters, sex) %>%
  table() %>%
  prop.table(margin = 1) %>%
  .[ubc_subclusters, ]
dset <- human_ubcs[[]] %>%
  select(subclusters, dataset_name) %>%
  mutate(dataset_name = str_remove(dataset_name, "_full_cerebellum_human")) %>%
  table() %>%
  prop.table(margin = 1) %>%
  .[ubc_subclusters, ]
pal1 <- RColorBrewer::brewer.pal(ncol(ages), "YlOrRd") # age
pal2 <- c("pink","blue") # sex
pal3 <- RColorBrewer::brewer.pal(3, "Dark2") # dataset

topannot <- HeatmapAnnotation(
  age = anno_barplot(
    ages, 
    gp = gpar(fill = pal1, border="white", lty="solid"),
    bar_width = 1, 
    height = unit(2, "in"),
    axis=FALSE
  ),
  sex = anno_barplot(
    sex,
    gp = gpar(fill = pal2, border="white"),
    bar_width = 1,
    height = unit(1,"in"),
    axis=FALSE
  ),
  dataset = anno_barplot(
    dset,
    gp = gpar(fill = pal3, border="white"),
    bar_width = 1,
    height = unit(1,"in"),
    axis=FALSE
  )
)
lgd_list <- list(
  Legend(labels=colnames(ages),title="age",legend_gp=gpar(fill=pal1)),
  Legend(labels=colnames(sex),title="sex",legend_gp=gpar(fill=pal2)),
  Legend(labels=colnames(dset),title="dataset",legend_gp=gpar(fill=pal3))
)

# genes to annotate
annoGenes <- map(
  .x = ranked_markers,
  .f = \(df) {
    head(df$feature, 5)
  }
) %>%
  unlist() %>%
  unname() %>%
  c("OTX2", "CBFA2T2", "SOX4", "SOX11", "LMX1A", "EOMES")


## -----------------------------------------------------------------------------
# get data slot, then manually scale the results
xpr <- AverageExpression(
  human_ubcs,
  assays = "SCT",
  features = unique_markers,
  group.by = "subclusters",
  slot = "data"
)[[1]]
xpr <- t(scale(t(xpr)))
# change column order
colnames(xpr) <- sub("-","_", colnames(xpr))
xpr <- xpr[, ubc_subclusters]

# heatmap colours
col_fun <- circlize::colorRamp2(
  breaks = c(-2,0,2),
  colors = c("#FF00FF","#000000", "#FFFF00")
)

ha <- rowAnnotation(
  foo = anno_mark(at = match(annoGenes, rownames(xpr)), labels = annoGenes)
)

hm <- Heatmap(
  xpr,
  name="scaled\nexpression",
  cluster_rows=TRUE,
  cluster_columns=FALSE,
  col = col_fun,
  show_row_names = FALSE,
  column_title="SCT data slot then manual scale()",
  top_annotation = topannot,
  right_annotation = ha
)

pdf(
  file= file.path(outDir, sprintf("UBCheatmap_top%d_genes.pdf", topX)),
  width = 8,
  height = 12
)
draw(hm, annotation_legend_list = lgd_list, merge_legend = TRUE)
dev.off()

