# plot UBC proportion 
rm(list=ls())

library(ggplot2)
library(tidytable)
library(readr)
library(stringr)
library(tibble)
library(Seurat)
library(harmony)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(BiocParallel)
library(ggrepel)


thesis_dir <- "../../"
source(file.path(thesis_dir, "./utils.R"))
source(file.path(thesis_dir, "../software/utilities/plotting.R"))

rpcaFile <- "/home/rstudio/isilon/private/llau/results/integrated/20240613/rpca/combined_clusters_rl.qs"
out_dir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/RPCA"
my_pals <- get_custom_pals()

dt <- format(Sys.Date(),"%y%m%d")
out_dir <- sprintf("%s/%s", out_dir, dt)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = FALSE)
}

logFile <- sprintf("%s/%s_rpca_ubc_clusters.log", out_dir, dt)
sink(logFile,split=TRUE)

tryCatch({
cat("Reading RPCA file\n")
t0 <- Sys.time()
rpca <- qs::qread(rpcaFile)
print(Sys.time() - t0)

cat("Plotting clusters of full rhombic lip RPCA integration\n")
p <- DimPlot(rpca, group.by="new_clusters")
ggsave(filename = "rpca_ubc_clusters.pdf", plot = p, path = out_dir, width = 6, height = 5)

#### do a dot plot of rhombic lip markers
###RLmarkers <- c("WLS","MKI67","SOX2","EOMES","PAX6",
###  "RBFOX3","ITPR1","CA8","AQP4","GFAP","OTX2","SOX4")
###p <- DotPlot(rpca, features = RLmarkers, group.by = "new_clusters") 
###ggsave(
###  filename = "rpca_rlmarkers_dotplot.pdf",
###  plot = p, path = out_dir, width = 10,height = 5)
###
cat("\n\n********************************\n")
cat(" Subsetting for UBC clusters\n")
ubc_clusters <- c("UBC_0", "UBC_1", "UBC_2", "UBC_3", "UBC_4", "UBC_5","UBC_6")
rpca_ubc <- subset(rpca, new_clusters %in% ubc_clusters)
cat(sprintf("Number of UBC cells: %d\n", ncol(rpca_ubc)))

cat("\n\n********************************\n")
cat(" Finding clusters with UBC subset\n")
DefaultAssay(rpca_ubc) <- "integrated"
rpca_ubc <- RunPCA(rpca_ubc)
rpca_ubc <- RunUMAP(rpca_ubc, dims = 1:30)
rpca_ubc <- FindNeighbors(rpca_ubc, dims = 1:30)
rpca_ubc <- FindClusters(rpca_ubc, resolution = 0.2)
p <- DimPlot(rpca_ubc, group.by = "integrated_snn_res.0.2") + theme_minimal(base_size = 18)
ggsave(filename = "rpca_ubc_integrated.pdf", plot = p, path = out_dir, width = 6, height = 5)
rpca_ubc$seurat_clusters <- rpca_ubc$integrated_snn_res.0.2

browser()

cat("\n\n********************************\n")
cat(" Creating Milo object\n")
srat <- rpca_ubc
md <- data.frame(srat[[]])
x <- mutate(
  md,
  sample = case_when(
    str_detect(dataset_name, "^Aldinger") ~ sample_id,
    .default = orig.ident
  ),
  cluster = new_clusters
) %>%
  select(cluster, sample, species, dataset_name)

# create a species column such that mouse becomes 1_mouse and human becomes 2_human
srat$species_num <- case_when(
  srat$species == "mouse" ~ "1_mouse",
  srat$species == "human" ~ "2_human",
  TRUE ~ as.character(srat$species)
)
srat$cluster <- x$cluster
srat$sample <- x$sample

cat("converting to SingleCellExperiment\n")
ubc_sc <- as.SingleCellExperiment(srat)
cat("creating milo object\n")
ubc_milo <- Milo(ubc_sc)

cat("\n\n********************************\n")
cat(" Building graph and testing for differential abundance with Milo\n")
d <- 30
mc <- MulticoreParam(workers = 8)
for (k in 50) {
  cat("********\n")
  cat("processing k=", k, "\n")
  cat("********\n")

  kDir <- sprintf("%s/k%d_d%d", out_dir, k, d)
  if (!dir.exists(kDir)) {
  dir.create(kDir, recursive = FALSE)
  }

  # took 17 minutes with 8 workers
  cat(sprintf("Building graph with k=%d, d=%d\n", k, d))
  t0 <- Sys.time()
  ubc_milo <- buildGraph(ubc_milo, k = k, d = d) #, BPPARAM = mc)
  print(Sys.time() - t0)

  cat("Calculating neighborhood abundance\n")
  t0 <- Sys.time()
  ubc_milo <- makeNhoods(ubc_milo, prop = 0.1, k = k, d=d, refined = TRUE)
  print(Sys.time() - t0)

  # print number of neighbourhoods
  cat("number of neighbourhoods:", nrow(ubc_milo@nhoods), "\n")

  cat("plotting neighborhood size histogram\n")
  p <- plotNhoodSizeHist(ubc_milo) + theme_minimal(base_size = 18) 
  ggsave(sprintf("%s/nbhd_size_hist.pdf", kDir), p, width = 4, height = 4)

  cat("counting cells\n")
  ubc_milo <- countCells(ubc_milo, meta.data = x, sample="sample")

  traj_design <- data.frame(colData(ubc_milo))[,c("sample", "species_num")]
  traj_design <- distinct(traj_design)
  traj_design
  rownames(traj_design) <- traj_design$sample

  cat("Calculating neighborhood distances\n")
  t0 <- Sys.time()
  ubc_milo <- calcNhoodDistance(ubc_milo, d=d)
  print(Sys.time() - t0)
  
  cat("Testing neighborhood differential abundance\n")
  da_results <- testNhoods(ubc_milo, 
    design = ~ species_num, 
    design.df = traj_design
  )

  write.table(da_results, file = sprintf("%s/da_results.txt", kDir), 
    sep = "\t", quote = FALSE, row.names = TRUE)
  ubc_milo <- buildNhoodGraph(ubc_milo)  

  cat("Plotting\n")
   lv <- levels(ubc_milo$seurat_clusters)
   pal <- setNames(
        pals::brewer.set2(length(lv)),
        nm = lv #paste0("UBC_", c(0:6))
  )
  p <- plotUMAP(ubc_milo, colour_by = "seurat_clusters")
  p <- p + scale_colour_manual(values = pal) + 
    theme_minimal(base_size = 18)
  p <- p + plotNhoodGraphDA(ubc_milo, da_results, alpha=0.05) 
  p <- p + theme(text = element_text(size = 18))
  ggsave(sprintf("%s/nbhd_umap.pdf", kDir), p, width = 15, height = 6)

}

cat("\n\n********************************\n")
cat(" Plotting cluster proportions\n")
# plot proportion of datasets in each cluster
md <- srat[[]]
  p <- ggplot(md, aes(x=seurat_clusters, fill=dataset_name)) +
    geom_bar(position = "fill") +
    theme_minimal(base_size = 24) +
    labs(y = "Proportion of cells", x = "Cluster", fill = "Dataset") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf("%s/cluster_datasets.pdf", out_dir), p, width = 10, height = 6)

p <- ggplot(md, aes(x=seurat_clusters, fill=species)) +
    geom_bar(position = "fill") +
    theme_minimal(base_size = 24) +
    labs(y = "Proportion of cells", x = "Cluster", fill = "Species") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf("%s/cluster_species.pdf", out_dir), p, width = 10, height = 6)

cat("\n\n********************************\n")
cat(" Run DEG analysis to identify upregulated markers of UBC clusters with sufficient human cells\n")
# Run DEG analysis to identify upregulated markers of UBC clusters with sufficient
# human cells 
human <- subset(srat, species == "human")
cat(sprintf("Number of human cells: %d\n", ncol(human)))

human.split <- SplitObject(human, split.by = "dataset_name")
human.split <- lapply(human.split, function(x) {
  x <- SCTransform(x, vst.flavor="v2", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)
})

cat("\nFinding integration features and anchors for human UBC clusters\n")
features <- SelectIntegrationFeatures(human.split, nfeatures = 3000)
human.split <- PrepSCTIntegration(human.split, anchor.features = features)
anchors <- FindIntegrationAnchors(human.split, normalization.method = "SCT", anchor.features = features)
human <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
cat(sprintf("Number of human cells: %d\n", ncol(human)))
# run DEG analysis on human UBC clusters
# Run PrepSCTFindMarkers
DefaultAssay(human) <- "SCT"
human <- PrepSCTFindMarkers(human, assay = "SCT")
source("clusterUtils.R")
Idents(human) <- "seurat_clusters"
deg <- run_all_de(
  srat = human,
  clusters = c(0:2,4),
  grouping.var = "dataset_name"
)
write.table(deg, file = sprintf("%s/human_ubc_deg.tsv", out_dir), sep = "\t", quote = FALSE, row.names = TRUE)

# plot a volcano plot for each cluster tested. label the top genes in each.
# color upregulated genes in red and downregulated genes in blue. use a significance threshold of 0.05 for adjusted p-value and a log2 fold change threshold of 0.25 for labeling genes.
for (cl in unique(deg$seurat_cluster)) {
    deg_cl <- deg %>% filter(seurat_cluster == cl)
     deg_cl <- deg_cl %>%
      mutate(
        significant = case_when(
          p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "upregulated",
          p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "downregulated",
          TRUE ~ "not_significant"
        )
      )
        p <- ggplot(deg_cl, aes(x=avg_log2FC, y=-log10(p_val_adj), color=significant)) +
      geom_point(size=0.5) +
      scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "not_significant" = "grey")) +
      labs(x = "Average log2 fold change", y = "-log10 adjusted p-value") +
      ggtitle(sprintf("Cluster %s", cl)) +
      theme_minimal(base_size = 18) +
      theme(legend.position = "none")
    # label top 50 genes by adjusted p-value
    top_genes <- deg_cl %>% arrange(p_val_adj) %>% head(50) %>% pull(gene)
    # use ggrepel to label the top genes on the volcano plot
    p <- p + geom_text_repel(data = deg_cl %>% filter(gene %in% top_genes), 
      aes(label = gene), size = 3, max.overlaps = Inf)
    ggsave(sprintf("%s/volcano_human_cluster_%s.pdf", out_dir, cl), plot = p, width = 6, height = 5)
}

# plot Violin Plots of EOMES for rpca_ubc
DefaultAssay(rpca_ubc) <- "RNA"
p <- VlnPlot(rpca_ubc, features = "EOMES", group.by = "seurat_clusters") + 
  theme_minimal(base_size = 18) +
  labs(x = "Cluster", y = "EOMES expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf("%s/EOMES_violin.pdf", out_dir), plot = p, width = 6, height = 5)


}, error=function(ex){
  print(ex)
}, finally={
  sink(NULL)
})