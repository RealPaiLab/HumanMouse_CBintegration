# Plot QC metrics for human and mouse datasets used for inference of UBC clusters
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(Seurat)
library(ggplot2)

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/QCmetrics"
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE);    
}
inFile <- "/home/rstudio/isilon/private/llau/results/integrated/20240524/25_pc_without_luo/25_pc_rl.qs"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = FALSE);
}


srat <- qs::qread(inFile)
srat$percent.mt <- PercentageFeatureSet(srat, pattern = "^MT-")

#' plot violin plots of percent.mt grouped by dataset_name
md <- srat[[]]
md$dataset_name <- sub("_full_cerebellum_", " ", md$dataset_name)

# add cell count and skinny boxplot to each violin
# remove labels beneath each violin plot
p <- ggplot(md, aes(x = dataset_name, y = percent.mt, fill = dataset_name)) +
  geom_violin() + geom_boxplot(width=0.1, cex=0.5) + 
  theme_minimal(base_size = 20) + theme(axis.text.x=element_blank())+
  labs(x = "Dataset", y = "% Mitochondrial genes")
ggsave(filename = sprintf("%s/percent_mt_violin.pdf", outDir), plot = p, width = 8, height = 4)

# now do the same for nFeature_RNA
p <- ggplot(md, aes(x = dataset_name, y = nFeature_RNA, fill = dataset_name)) +
  geom_violin() + geom_boxplot(width=0.1, cex=0.5) +
  theme_minimal(base_size = 20) + theme(axis.text.x=element_blank()) + 
  labs(x = "Dataset", y = "Number of Features (Genes)")
ggsave(filename = sprintf("%s/nFeature_RNA_violin.pdf", outDir), plot = p, width = 8, height = 4)

# now do the same for nCount_RNA
p <- ggplot(md, aes(x = dataset_name, y = nCount_RNA, fill = dataset_name)) +
  geom_violin() + geom_boxplot(width=0.1, cex=0.5) +
  theme_minimal(base_size = 20) + theme(axis.text.x=element_blank()) + 
  labs(x = "Dataset", y = "Number of Molecules")
ggsave(filename = sprintf("%s/nCount_RNA_violin.pdf", outDir), plot = p, width = 8, height = 4)

