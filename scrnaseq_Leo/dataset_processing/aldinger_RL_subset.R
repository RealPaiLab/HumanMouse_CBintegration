library(Seurat)
library(tidyverse)
library(qs)

setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")
source("scrnaseq_Leo/utilities/isolate_RL_lineage.R")

# Luo human FC
isolate_RL_lineage("/.mounts/labs/pailab/private/llau/data/Aldinger_2021/Aldinger_full_cerebellum_human_20240325_stand.qs", 
                    "/.mounts/labs/pailab/private/llau/data/Aldinger_2021", 
                    "Aldinger", 
                    "human")