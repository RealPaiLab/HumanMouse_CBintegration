
# Milo testing of human-mouse UBCs
# following tutorial from : https://marionilab.github.io/miloR/articles/milo_demo.html#1-defining-representative-neighbourhoods

rm(list=ls())

library(miloR)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(Seurat)
library(stringr)
library(tidyr)
library(scater)

cat("reading qs\n")
t0 <- Sys.time()
srat <- qs::qread("/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20240825/ubc_subset.qs")
print(Sys.time() - t0)

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/Milo"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

logFile <- sprintf("%s/milo_log.txt", outDir)
sink(logFile, split=TRUE)

tryCatch({


md <- data.frame(srat[[]])
x <- mutate(
  md,
  sample = case_when(
    str_detect(dataset_name, "^Aldinger") ~ sample_id,
    .default = orig.ident
  ),
  cluster = snn_res.0.3
) %>%
  select(cluster, sample, species, dataset_name)

srat$sample <- x$sample
srat$cluster <- paste0("UBC_", x$cluster)

# create a species column such that mouse becomes 1_mouse and human becomes 2_human
srat$species_num <- case_when(
  srat$species == "mouse" ~ "1_mouse",
  srat$species == "human" ~ "2_human",
  TRUE ~ as.character(srat$species)
)

#### run SCT, PCA, and plot UMAP on the subsetted data
###srat <- SCTransform(srat, verbose = FALSE)
###srat <- RunPCA(srat, ndims = 30, verbose = FALSE)
###srat <- RunUMAP(srat, dims = 1:15, verbose = FALSE)
#### plot UMAP colored by cluster
p <- DimPlot(srat, group.by = "cluster", label = TRUE) + 
  theme(legend.position = "none")
ggsave(sprintf("%s/seurat_umap_clusters.png", outDir), p, 
width = 6, height = 6)

cat("converting to SingleCellExperiment\n")
ubc_sc <- as.SingleCellExperiment(srat)
cat("creating milo object\n")
ubc_milo <- Milo(ubc_sc)

# k=10 didn't work, got an NA error in calcNhoodDistance.

for (k in c(40, 50, 70)) { #c(35, 40, 45, 50)) { #15,25,30)) {
cat("********\n")
cat("processing k=", k, "\n")
cat("********\n")
# k of 30 selected to ensure neighbourhood sizes peak between 50 and 100 
# as suggested in this discussion from the miloR developers:
# https://github.com/MarioniLab/miloR/issues/350
d <- 30

kDir <- sprintf("%s/k%d_d%d", outDir, k, d)
if (!dir.exists(kDir)) {
  dir.create(kDir, recursive = FALSE)
}

cat(sprintf("Building graph with k=%d, d=%d\n", k, d))
ubc_milo <- buildGraph(ubc_milo, k = k, d = d)
cat("Calculating neighborhood abundance\n")
ubc_milo <- makeNhoods(ubc_milo, prop = 0.1, k = k, d=d, refined = TRUE)
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

# results for testNHoods interpretation
# https://github.com/MarioniLab/miloR/issues/81
# logFC is the log fold change between conditions. For a comparison across ordered variables,
### as is the case for E7 - E8, then the interpretation of the logFC is the mean per-unit 
### change in cell counts per change in predictor variable. 
### So in the example above the first line of the table shows that there is an 
### average 7-fold lower counts as you go from E7 to E7.5 to E8. In essence, 
### you interpret it in the same way that you would a regression coefficient from 
### any other linear model.
#
t0 <- Sys.time()
ubc_milo <- calcNhoodDistance(ubc_milo, d=d)
print(Sys.time() - t0)
da_results <- testNhoods(ubc_milo, 
  design = ~ species_num, 
  design.df = traj_design)

write.table(da_results, file = sprintf("%s/da_results.txt", kDir), 
    sep = "\t", quote = FALSE, row.names = TRUE)

pal <- setNames(
      pals::brewer.set2(6),
      nm = paste0("UBC_", c(0:5))
    )
ubc_milo <- buildNhoodGraph(ubc_milo)

cat("Plotting\n")
p <- plotUMAP(ubc_milo, colour_by = "cluster")
p <- p + scale_colour_manual(values = pal) + 
  theme_minimal(base_size = 18)
p <- p + plotNhoodGraphDA(ubc_milo, da_results, alpha=0.05) 
p <- p + theme(text = element_text(size = 18))
ggsave(sprintf("%s/nbhd_umap.pdf", kDir), p, width = 15, height = 6)

# get a list of significantly upregulated neighbourhoods
sig_nhoods <- da_results$Nhood[da_results$SpatialFDR < 0.05 & da_results$logFC > 0]
cat(sprintf("number of significant neighbourhoods: %d of %d\n", length(sig_nhoods), nrow(da_results)))

cell_member <- as.matrix(ubc_milo@nhoods)
nindex <- unlist(ubc_milo@nhoodIndex)
sig_nhoods_idx <- nindex[sig_nhoods]

colidx <- which(colnames(cell_member) %in% sig_nhoods_idx)
in_sig_nhoods <- rowSums(cell_member[, colidx, drop=FALSE])
cells_in_sig_nhoods <- rownames(cell_member)[in_sig_nhoods > 0]
cat(sprintf("number of cells in significant neighbourhoods: %d\n", 
  length(cells_in_sig_nhoods)))

# calculate, for each subcluster, what fraction of cells are in significant neighbourhoods, and plot as a barplot
md <- data.frame(colData(ubc_milo))
out <- list()
for (cl in unique(md$subclusters)) {
  cells_in_cl <- rownames(md)[md$subclusters == cl]
  n_cells_in_cl <- length(cells_in_cl)
  n_cells_in_sig_nhoods_in_cl <- sum(cells_in_cl %in% cells_in_sig_nhoods)
  prop_in_sig_nhoods <- n_cells_in_sig_nhoods_in_cl / n_cells_in_cl
  out[[cl]] <- data.frame(
    subcluster = cl,
    n_cells = n_cells_in_cl,
    n_cells_in_sig_nhoods = n_cells_in_sig_nhoods_in_cl,
    prop_in_sig_nhoods = prop_in_sig_nhoods
  )
}
out_df <- do.call(rbind, out)

# show as a percentage
p <- ggplot(out_df, aes(x = subcluster, y = prop_in_sig_nhoods * 100, fill = subcluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none") +
  labs(x = "Subcluster", y = "% cells in enriched neighbourhoods")
ggsave(sprintf("%s/prop_cells_in_sig_nhoods_by_subcluster.pdf", kDir), p, width = 6, height = 4)

}
cat("done\n")


}, error = function(e) {
  print(e)
}, finally = {
  sink()
})