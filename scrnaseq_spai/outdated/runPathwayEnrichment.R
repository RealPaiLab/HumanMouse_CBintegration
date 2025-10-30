rm(list=ls())
library(tidyverse)
library(gprofiler2)

# load helper functions
source("/home/rstudio/isilon/private/icheong/CBL_scRNAseq/software/utilities/gprofiler2_helpers.R")

inRoot <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/UBCclusters_241112"

for (curRes in c(0.6)) {
    DEGresults <- sprintf("%s/%1.1f_res/%1.1f_markers.csv",
        inRoot,curRes,curRes)

    outDir <- dirname(DEGresults)

    cat(sprintf("%1.2f: reading DEG\n",curRes))
    res <- read.delim(DEGresults,sep=",",h=TRUE)

    for (k in unique(res$cluster)) {
        tmp <- subset(res, cluster %in% k)
        bg <- unique(tmp$feature)
        fg <- unique(tmp$feature[tmp$p_val < 0.05])
        cat(sprintf("Cluster %i: %i fg, %i bg\n", k,length(fg),length(bg)))
    # run pathway enrichment
    gost_res <- run_gost(
          query = fg,
          organism = "gp__OUFl_gpIE_RyE", # 10-250 genes in each pathway
          significant = FALSE,
          custom_bg = bg,
          filename =  sprintf("%s/cluster%i_pathwayEnrichment",outDir,k)
        )

    }
}