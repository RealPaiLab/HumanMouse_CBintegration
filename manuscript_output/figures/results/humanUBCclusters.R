rm(list=ls())
library(Seurat)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

srat_file <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBCclusters/fromQuang/UBC.Harmony.RDS"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/DimPlots"

DEdir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/diffExpr"
useDEfile <- sprintf("%s/250424/de_genes.tsv", DEdir)

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s",outDir,dt)

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

cat("Reading human UBC cluster file\n")
t0 <- Sys.time()
srat <- readRDS(srat_file)
print(Sys.time() - t0)

# remove cells with age "9 PCW" and "10 PCW"
srat <- subset(srat, subset = age != "9 PCW" & age != "10 PCW")

p <- DimPlot(srat, reduction = "umap", 
    group.by = "SCT_snn_res.0.5", 
    label = TRUE, label.size = 10, pt.size = 1.2,
    label.box=TRUE, repel = TRUE) +
    #cols=clrs$Colour[order(clrs$Cluster)]) +
  ggtitle(sprintf("Human UBC clusters, %i cells",ncol(srat))) 

p <- p + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))
ggsave(filename = sprintf("%s/UBCclusters_SCT_snn_res.0.5.pdf", outDir), 
    plot = p, width = 8, height = 6, dpi = 300)

# show a barplot of dataset_name by SCT_snn_res.0.5
p2 <- srat[[]] %>%
  group_by(SCT_snn_res.0.5, dataset_name) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = SCT_snn_res.0.5, y = count, fill = dataset_name)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) + 
  xlab("SCT_snn_res.0.5") + ylab("Cell count") +
  ggtitle("Dataset composition of human UBC clusters") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(filename = sprintf("%s/UBCclusters_SCT_snn_res.0.5_barplot.pdf", outDir), 
      plot = p2, width = 8, height = 6, dpi = 300)    