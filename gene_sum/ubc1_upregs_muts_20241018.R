# generate mutation counts for human upregulated UBC1 genes (comparing to all rest UBCs)
# definition provided by Shraddha in 20241017 2:14pm email
# genes that have (Sepp_full_cerebellum_human_p_val_adj < 0.05) && (Aldinger_full_cerebellum_human_p_val_adj < 0.05) && (Aldinger_full_cerebellum_human_avg_log2FC > 0) && (Sepp_full_cerebellum_human_avg_log2FC > 0)
# also forced genes that might be of interests: c("CNTNAP2", "ARHGAP11B", "FOXP2", "RRM2", "EOMES")


rm(list=ls())

library(GenomicRanges)
library(rtracklayer)
library(dplyr)

default_wd <- this.path::this.dir()
setwd(default_wd) # set current scripts' dir as working dir
source("./plot_geneExploration.R")
source("./utils.R")

### CONFIGS ###
## INPUT ##
# GENCODE #
gencode_file <- "/.mounts/labs/pailab/private/projects/FetalHindbrain/anno/gencode.v44.basic.annotation.gtf"

# G3&4 MB PCAWG SNV/INDEL in hg38 #
# Using PCAWG (downloaded from ICGC portal) because the high calling quality
# Also can change to PEMECA, PEMECA-PCAWG if needed
mb_mut_file <- list(
  "snv-indel" = list(
    g3_mut_file = "/data/xsun/20240314/mutation/PCAWG/merged/Group3_PCAWG_snv-indel_hg38.bed",
    g4_mut_file = "/data/xsun/20240314/mutation/PCAWG/merged/Group4_PCAWG_snv-indel_hg38.bed",
    shhmut_file = "/data/xsun/20240314/mutation/PCAWG/merged/SHH_PCAWG_snv-indel_hg38.bed"
  ),
  "snv-indel-sv" = list(
    g3_mut_file = "/data/xsun/20241005_PCAWG_MB_muts/Group3_PCAWG_snv-indel-sv_hg38.bed",
    g4_mut_file = "/data/xsun/20241005_PCAWG_MB_muts/Group4_PCAWG_snv-indel-sv_hg38.bed",
    shhmut_file = "/data/xsun/20241005_PCAWG_MB_muts/SHH_PCAWG_snv-indel-sv_hg38.bed"
  )
)

## OUTPUT ##
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf(
  #"/.mounts/labs/pailab/private/projects/HumanMouseUBC/TumourAnalysis/ubc1_upregs_muts/%s",
  "/.mounts/labs/pailab/private/xsun/output/HumanMouseUBC/TumourAnalysis/ubc1_upregs_muts/%s",
  dt
  )

if (! dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}


### MAIN ###
## PREP ##
# gencode #
gencode_gff <- rtracklayer::readGFF(gencode_file)
gencode <- GenomicRanges::makeGRangesFromDataFrame(gencode_gff, keep.extra.columns = T)

# mutation #
# with sv
mut_class <- "snv-indel-sv"

g3_mut <- import.bed(mb_mut_file[[mut_class]]$g3_mut_file)
g3_mut$class <- "Grp3"

g4_mut <- import.bed(mb_mut_file[[mut_class]]$g4_mut_file)
g4_mut$class <- "Grp4"

shhmut <- import.bed(mb_mut_file[[mut_class]]$shhmut_file)
shhmut$class <- "SHH"

combined_mut <- c(g3_mut, g4_mut, shhmut)

# without sv
mut_class <- "snv-indel"

g3_mut <- import.bed(mb_mut_file[[mut_class]]$g3_mut_file)
g3_mut$class <- "Grp3"

g4_mut <- import.bed(mb_mut_file[[mut_class]]$g4_mut_file)
g4_mut$class <- "Grp4"

shhmut <- import.bed(mb_mut_file[[mut_class]]$shhmut_file)
shhmut$class <- "SHH"

combined_mut_snvindel <- c(g3_mut, g4_mut, shhmut)


# count muts for all UBC1 up TFs #
# copied from 
# /.mounts/labs/pailab/private/llau/results/integrated/20240715/all_tested_genes.csv
# /.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240827/diff_expr_markers/cluster_marker_ranking.rds

ubc1_root <- "/data/xsun/ubc1"

# load differential gene expression results (all UBC clusters, both datasets)
ubc_diff_expr <- read.csv(sprintf("%s/all_tested_genes.csv", ubc1_root)) %>%
  rename_with(.fn = stringr::str_remove, pattern = "full_cerebellum_")

# select genes
addional_genes <- c("CNTNAP2", "ARHGAP11B", "FOXP2", "RRM2", "EOMES")

interested_genes <- ubc_diff_expr %>%
  filter(cluster == "UBC_1") %>%
  filter(Sepp_human_p_val_adj < .05 &
           Sepp_human_avg_log2FC > 0 &
           Aldinger_human_p_val_adj < .05 &
           Aldinger_human_avg_log2FC > 0) %>%
  pull(feature)


# Count muts #
cluster <- "UBC1"
dict <- list(UBC1 = unique(c(interested_genes, addional_genes)))

# generate dfs
# snv/indel/sv
snvindelsv_df <- generate_mut_df(dict$UBC1, 
                                 mutation = combined_mut, 
                                 gencode = gencode)
if (sum(is.na(res)) > 0) {
  warning("--- There are genes not found in Gencode! Remove/rename the genes!")
}

write.table(snvindelsv_df, 
            sprintf("%s/UBC1_upregs_muts_snv-indel-sv.tsv", outDir), 
            col.names = T, row.names = T, quote = F, sep = "\t"
            )

# snv/indel
snvindel_df <- generate_mut_df(dict$UBC1, 
                               mutation = combined_mut_snvindel, 
                               gencode = gencode)

write.table(snvindel_df, 
            sprintf("%s/UBC1_upregs_muts_snv-indel.tsv", outDir), 
            col.names = T, row.names = T, quote = F, sep = "\t"
)

# pyramid plot
.plot <- plot_mut_pyramid(cluster = "UBC1", 
                 dict = dict, 
                 combined_mut_snvindel = combined_mut_snvindel, 
                 combined_mut_snvindelsv = combined_mut, 
                 gencode = gencode
                 )

ggsave(sprintf("%s/UBC1_upregs_muts.png", outDir), .plot, 
       width = 23, height = 18, dpi = 600)


























