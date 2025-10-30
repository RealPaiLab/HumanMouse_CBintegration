library(Seurat)
library(tidyverse)
library(qs)

setwd("/u/llau/software/mb_scrnaseq/MB_scRNAseq")
source("scrnaseq_Leo/utilities/isolate_RL_lineage.R")

# Sepp human FC
isolate_RL_lineage("/.mounts/labs/pailab/private/llau/data/Sepp_2024/Sepp_full_cerebellum_human_20240404_stand.qs", 
                    "/.mounts/labs/pailab/private/llau/data/Sepp_2024", 
                    "Sepp", 
                    "human")

# Sepp mouse FC
isolate_RL_lineage("/.mounts/labs/pailab/private/llau/data/Sepp_2024/Sepp_full_cerebellum_mouse_20240404_stand.qs", 
                    "/.mounts/labs/pailab/private/llau/data/Sepp_2024", 
                    "Sepp", 
                    "mouse")
