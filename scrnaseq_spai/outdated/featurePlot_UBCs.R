rm(list=ls())

library(tidyverse)
library(patchwork)
library(Seurat)
library(pals)
library(readxl)

srat_qs <- get_srat_paths()
ubc_srat <- load_srat(srat_qs["ubc"]) %>%
  pluck("ubc")
##message("done")
my_pals <- get_custom_pals()


humanUBC <- subset(ubc_srat, subset = species == "human")
DefaultAssay(humanUBC) <- "RNA"

humanUBC <- NormalizeData(humanUBC) 
humanUBC <- SCTransform(humanUBC, vars.to.regress="CC.Difference") 
# UMAP viz & DEG
0

browser()
# https://satijalab.org/seurat/articles/essential_commands.html
# head(humanUBC[[]])
# head(humanUBC@meta.data)
# head(humanUBC$subclusters)
# get normalized expressed: xpr <- humanUBC[["RNA"]]$data

# ASSAYS
# 1. RNA - counts: raw counts, data: normalized counts, scaled.data: normalized and scaled to (mean 0, SD 1)
# 2. SCT - results of SCTransform - used for viz of UMAP.

# Get per-cluster expression of FOXP2
# note that this does an mean(expm1(RNA counts))
# it assumes LogNormalize has been run on the RNA data.
# avg <- AverageExpression(humanUBC, 
#   features="FOXP2", 
#   group.by="subclusters",assays=)

# todo - try assay=RNA and assay=SCT assay
###avg <- AverageExpression(humanUBC, 
###  features=huttner57, 
###  group.by="subclusters",
###  assays = "RNA"
###  )

# 1. todo - try RNA and SCT for dot plot.
# what is DotPlot cutoff for "expressed" in a cell? Is it any > 0?
DefaultAssay(humanUBC) <- "RNA"
for (assayType in c("RNA","SCT")) {
  p <- DotPlot(humanUBC, 
    features=c("FOXP2","OTX2","EOMES",huttner57[1:10]),
    group.by = "subclusters",
    assay = assayType) # try SCT as well. 
  ggsave(p,file=sprintf("humanUBC_dotPlot_%s.png",assayType),
    path=out_dir, 
    width=13, height=6,units="in")
}


subcluster_umap <- DimPlot(
  ubc_srat,
  reduction = "umap",
  group.by = "subclusters",
  cols = my_pals$ubc_integ_clust,
  label = TRUE,
  label.size = 3,
  repel = TRUE
) + 
  labs(title = NULL) + 
  theme_classic2() + 
  theme(legend.position = "none")

# hacky way of setting the seed for the UMAP labels
subcluster_umap[[1]]$layers[[2]]$geom_params$seed <- 42

##ggsave(
##  filename = "subcluster_umap.png",
##  plot = subcluster_umap,
##  path = out_dir,
##  width = 5,
##  height = 5,
##  units = "in",
##  dpi = 300
##)

geneSet <-  c("EOMES","FOXP2","CNTNAP2","OTX2", "RBFOX1")
p <- FeaturePlot(ubc_srat, 
    geneSet,
    order=TRUE,
    ncol=3)
ggsave(p, 
    file = "FeaturePlot.png",
    width = 11,
    height = 8,
    units = "in",
    path = out_dir
    )


# plot density plot for FeaturePlot
for (curGene in geneSet) {
    df <- get_expression_data(ubc_srat, genes = curGene)
    umap_cor<- ubc_srat@reductions$umap@cell.embeddings  %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "cell")

    df<- left_join(df, umap_cor)
    colnames(df)[2] <- "GENE"

    p3<- ggplot(df, aes(x= UMAP_1, y= UMAP_2 )) +
            geom_point(data = df %>% filter(GENE == 0), 
                color = "#440154FF", size = 0.6) +
            ggpointdensity::geom_pointdensity(
                data = df %>% filter(GENE > 0), size = 0.6) +
            viridis::scale_color_viridis() +
            theme_classic(base_size = 14) +
            ggtitle(curGene)
    ggsave(p3,
        file=sprintf("DensityPlot_%s.png",curGene),
        path=out_dir)
    rm(p3)
}

p <- DotPlot(ubc_srat, features=geneSet, split.by="subclusters")
ggsave(p,file="DotPlot.png",path=out_dir)