library(GenomicRanges)
library(rtracklayer)

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


# G3&4 MB bulk RNA-seq from Hendrikse et al #
hendrikse_meta_file <- "/.mounts/labs/pailab/src/medulloblastoma-genomics/RNA/Hendrikse_2022/41586_2022_5215_MOESM4_ESM.xlsx"
hendrikse_expr_file <- "/.mounts/labs/pailab/src/medulloblastoma-genomics/RNA/Hendrikse_2022/G3_G4_raw_counts_bulkRNA_Hendrikse2022.txt"

# G3&4 MB bulk RNA micro-array from Cavali et al #
cavalli_meta_file <- "/.mounts/labs/pailab/src/medulloblastoma-genomics/RNA/Cavalli_2017/mmc2.xlsx"
cavalli_expr_file <- "/.mounts/labs/pailab/src/medulloblastoma-genomics/RNA/Cavalli_2017/GSE85217_M_exp_763_MB_SubtypeStudy_TaylorLab.txt.gz"


## OUTPUT ##
outDir <- "/data/xsun/tmp"
dt <- format(Sys.Date(),"%y%m%d")

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

# cavalli RNA-seq #
# meta
cavalli_meta <- get_cavalli_meta(cavalli_meta_file)

# expr
cavalli_expr <- get_cavalli_expr(cavalli_expr_file)

# hendrikse RNA-seq #
hendrikse_dds <- get_hendrikse_dds(hendrikse_meta_file, hendrikse_expr_file)



## PLOT ##
# EXAMPLES #
# Plot mutations
plot_mutation("BRCA1", combined_mut, gencode)$plot # No coding mut
plot_mutation("GSE1", combined_mut, gencode)$plot # SHH coding mut
plot_mutation("SMARCA4", combined_mut, gencode)$plot # Grp3&4 coding mut
plot_mutation("KDM6A", combined_mut, gencode)$plot # Grp4 coding mut
plot_mutation("TERT", combined_mut, gencode)$plot # SHH promoter mut


# plot expression
# Cavalli
plot_bulkRNA_cavalli(gene = "OTX2", # high in Grp3/4
                     meta_df = cavalli_meta, expr = cavalli_expr)
plot_bulkRNA_cavalli(gene = "MYC", # high in Grp3
                     meta_df = cavalli_meta, expr = cavalli_expr)
plot_bulkRNA_cavalli(gene = "ATOH1", # high in SHH
                     meta_df = cavalli_meta, expr = cavalli_expr)
plot_bulkRNA_cavalli(gene = "EOMES", 
                     meta_df = cavalli_meta, expr = cavalli_expr)

# hendrikse
plot_bulkRNA_hendrikse(gene = "OTX2", 
                     dds = hendrikse_dds)
plot_bulkRNA_hendrikse(gene = "MYC", # high in Grp3
                       dds = hendrikse_dds)
plot_bulkRNA_hendrikse(gene = "ATOH1", # high in SHH
                     dds = hendrikse_dds)
plot_bulkRNA_hendrikse(gene = "EOMES", 
                       dds = hendrikse_dds)

# Plot survival based on Cavalli
plot_survival("MYC", cavalli_meta, cavalli_expr)
plot_survival("MYC", cavalli_meta, cavalli_expr, subgroups = c("Group3"))
plot_survival("MYC", cavalli_meta, cavalli_expr, subgroups = c("Group4"))


# TARGETS #
# UBC_1 genes
ubc1_genes <- c("SOX4", "SOX11", "FOXP2", "EBF3", "ARID3A")
pdf(sprintf("%s/ubc1_genes_sum.pdf", outDir), width = 14, height = 10)
for (gene in ubc1_genes) {
  p <- plot_gene_sum(gene = gene, 
                      mutation = combined_mut, gencode = gencode,
                      meta_df = cavalli_meta, expr = cavalli_expr,
                      dds = hendrikse_dds, subgroups = c("Group4"))
  print(p)
}
dev.off()


# UBC_2 genes
ubc2_genes <- c("ZNF804A", "ZNF385D", "ZNF385B", "HIVEP2", "RFX3", "ZXDB")
pdf(sprintf("%s/ubc2_genes_sum.pdf", outDir), width = 14, height = 10)
for (gene in ubc2_genes) {
  p <- plot_gene_sum(gene = gene, 
                mutation = combined_mut, gencode = gencode,
                meta_df = cavalli_meta, expr = cavalli_expr,
                dds = hendrikse_dds, subgroups = c("Group4"))
  print(p)
}
dev.off()


# controls
control_genes <- c("TERT", "GSE1", "GLI2",
                   "OXT2", "MYC", "MYCN", 
                   "PRDM6", "CBFA2T2", "CBFA2T3", "RUNX1T1", "GFI1", "GFI1B", 
                   "KDM6A", "KDM2B"
                   )
pdf(sprintf("%s/control_genes_sum.pdf", outDir), width = 14, height = 10)
for (gene in ubc2_genes) {
  p <- plot_gene_sum(gene = gene, 
                     mutation = combined_mut, gencode = gencode,
                     meta_df = cavalli_meta, expr = cavalli_expr,
                     dds = hendrikse_dds, subgroups = c("Group4"))
  print(p)
}
dev.off()


# count muts for all UBC1 up TFs #
# copied from 
# /.mounts/labs/pailab/private/llau/results/integrated/20240715/all_tested_genes.csv
# /.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/integrated/20240827/diff_expr_markers/cluster_marker_ranking.rds

ubc1_root <- "/data/xsun/ubc1"

# load differential gene expression results (all UBC clusters, both datasets)
ubc_diff_expr <- read.csv(sprintf("%s/all_tested_genes.csv", ubc1_root)) %>%
  rename_with(.fn = stringr::str_remove, pattern = "full_cerebellum_")

# get significantly upregulated UBC markers (all UBC clusters, both datasets)
cluster_marker_ranking <- readRDS(sprintf("%s/cluster_marker_ranking.rds", ubc1_root))

# get the UBC clusters
ubc_subclusters <- names(cluster_marker_ranking)

# load list of TFs
tfs <- read.csv(
  "/.mounts/labs/pailab/src/transcription-factors/human-tfs-lambert/full_database_v1.01.csv",
  row.names = 1,
  check.names = FALSE
)

# get only confirmed TFs
tfs <- dplyr::filter(tfs, `Is TF?` == "Yes") %>%
  pull(`HGNC symbol`)
ubc_diff_expr <- dplyr::filter(
  ubc_diff_expr,
  feature %in% tfs
)

# get upregulated TFs
cluster_marker_ranking <- purrr::map(
  .x = cluster_marker_ranking,
  .f = \(x) {
    x <- dplyr::filter(x, feature %in% tfs)
    return(x)
  }
)

ubc1_up_tf <- cluster_marker_ranking$UBC_1$feature
ubc2_up_tf <- cluster_marker_ranking$UBC_2$feature


cluster <- "UBC1"
dict <- list(UBC1 = unique(c(ubc1_up_tf, 
                             c("CNTNAP2", "ARHGAP11B", "FOXP2", "RRM2", "EOMES")
                             )
                           ), 
             UBC2 = ubc2_up_tf)
plot_mut_pyramid("UBC1", dict, combined_mut_snvindel, combined_mut, gencode)

