
library(ggplot2)
library(ggrepel)
library(dplyr)

pySCENIC_tumour <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/pyscenic/20240928/vladoiu_mb_subtype_rss.csv"
pySCENIC_dev <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/pyscenic/20240704/aldinger_sepp_RL.regulons.dat"

deg_tumour <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/tumour/Vladoiu/20240917/mb_subtype_diff_exp.csv"
deg_rl <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/RLDEG/260309/RL_DEG.txt"


outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/CBintegrationPaper/RLDEG"
if (!file.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}
dt <- format(Sys.Date(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = FALSE)
}

logFile <- sprintf("%s/RL_to_Group3MB_analysis.txt", outDir)
sink(logFile, split=TRUE)
tryCatch({

# Plot regulons for Group 3 MB.
tum_reg <- read.delim(pySCENIC_tumour, sep = ",", header = TRUE)

g3 <- tum_reg[1,-1]; 
g3 <- as.numeric(g3); names(g3) <- colnames(tum_reg)[-1]
g3 <- g3[order(g3, decreasing = TRUE)]

g3 <- data.frame(regulon=names(g3), rank=1:length(g3), rss=g3)

# label the top 10 regulons
# colour the top 10 dots in the plot as well. when labelling the datapoints, make the text fan out to the right.
g3$label <- ifelse(g3$rank <= 10, g3$regulon, "")

clr <- "#56b4e9"
p <- ggplot(g3, aes(x=rank, y=rss)) +
  geom_point() +
  theme_classic() + theme(text = element_text(size=18), axis.text = element_text(size=24)) +
  labs(x="Regulon rank", y="Regulon specificity score (RSS)", 
  title="PySCENIC regulon ranking for Group 3 MB")
p <- p + geom_text_repel(aes(label=label), size=6, max.overlaps = Inf, col=clr, nudge_x = 10, direction = "y", hjust=0)
p <- p + geom_point(data=g3[g3$rank <= 10,], aes(x=rank, y=rss), color=clr, size=3)
ggsave(p,file=sprintf("%s/Group3MB_PySCENIC_regulon_ranking.pdf", outDir), width=6, height=6)

# Compare intersecting TFs upregulated in RL and Group 3 MB.
lambert_tf <- read.csv("/home/rstudio/isilon/src/transcription-factors/human-tfs-lambert/full_database_v1.01.csv") %>%
    filter(Is.TF. == "Yes") %>%
    select(tf = HGNC.symbol) %>%
    unique()
cat(sprintf("Read %i TFs from Lambert et al. human TF database\n", nrow(lambert_tf)))

deg_rl <- read.delim(deg_rl, sep="\t", header=TRUE)
deg_rl <- subset(deg_rl, seurat_cluster == "RL" & p_val_adj < 0.05 & avg_log2FC > 0)
cat(sprintf("%d genes upregulated in RL\n", nrow(deg_rl)))
deg_rl <- subset(deg_rl, gene %in% lambert_tf$tf)
cat(sprintf("%d of the genes upregulated in RL are TFs according to Lambert et al\n", nrow(deg_rl)))

deg_tumour <- read.delim(deg_tumour, sep=",", header=TRUE)
deg_tumour <- subset(deg_tumour, cluster == "G3"); deg_g3 <- deg_tumour
deg_tumour <- subset(deg_tumour, p_val_adj < 0.05 & avg_log2FC > 0)
cat(sprintf("%d genes upregulated in Group 3 MB\n", nrow(deg_tumour)))
deg_tumour <- subset(deg_tumour, gene %in% lambert_tf$tf)
cat(sprintf("%d of the genes upregulated in Group 3 MB are TFs according to Lambert et al\n", nrow(deg_tumour)))

# plot a Venn diagram of the overlap between the TFs upregulated in RL and Group 3 MB. label the TFs in the intersection.
library(VennDiagram)
tum_tfs <- deg_tumour$gene
rl_tfs <- deg_rl$gene
venn.plot <- draw.pairwise.venn(area1 = length(rl_tfs), area2 = length(tum_tfs), 
    cross.area = length(intersect(rl_tfs, tum_tfs)), 
    category = c("RL","Group 3 MB"), 
    fill = c( "#42d4f4","#56b4e9"), alpha = 0.5, cat.pos = c(0, 180), 
    cat.dist = 0.05, scaled = FALSE)
# plot and save as pdf
ggsave(venn.plot, file=sprintf("%s/RL_Group3MB_upregulated_TF_overlap.pdf", outDir), width=6, height=6)
# print the common TFs in the console
common_tfs <- intersect(rl_tfs, tum_tfs)
cat(sprintf("TFs upregulated in both RL and Group 3 MB: %s\n", paste(sort(common_tfs), collapse=", ")))

# plot volcano of Group 3 MB tumours.
source("clusterUtils.R")
deg_g3$seurat_cluster <- deg_g3$cluster
p <- plotVolcano_usingDEG(deg_g3, showTopGenes=15, logFCcutoff=0.25, titlePfx="Group 3 MB: ")
ggsave(p, file=sprintf("%s/Group3MB_volcano.pdf", outDir), width=10, height=6)

}, error=function(ex) {
 print(ex)
}, finally={
  sink()
})





