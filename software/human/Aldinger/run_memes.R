# ==============================================================================
# Motif enrichment testing with the MEME suite
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(argparse)
library(tidyverse)
library(GenomicRanges)
library(IRanges)
library(memes)
library(BSgenome.Hsapiens.UCSC.hg38)

# options(meme_path = "/home/rstudio/")

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # CSV of differentially expressed genes from `FindMarkers`
  "--deg_csv",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # use genes with logFC > 0
  "--logfc_pos",
  default = TRUE,
  action = "store_true"
)
parser$add_argument(
  # use genes with logFC < 0
  "--logfc_neg",
  dest = "logfc_pos",
  action = "store_false"
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # GENCODE GTF annotation file
  "--gencode_gtf",
  default = "/isilon/CBL_scRNAseq-archived/data/gencode/gencode.v42.basic.annotation.gtf"
)
parser$add_argument(
  # MEME Suite database
  "--meme_db",
  default = "/isilon/CBL_scRNAseq-archived/data/MEMES/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
)

if (interactive()) {
  args <- parser$parse_args(c(
    "--deg_csv", "/CBL_scRNAseq/results/human/Aldinger/20230205/de_genes.csv",
    "--logfc_pos",
    "--out_dir", "/CBL_scRNAseq/results/human/Aldinger/20230419/homol_ubc"
  ))
} else {
  args <- parser$parse_args()
}

if (is.null(args$out_dir)) {
  stop("Argument for `out_dir` is missing; please provide an output directory")
} else {
  out_dir <- args$out_dir
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
}

# ------------------------------------------------------------------------------
# import function

source("/CBL_scRNAseq/software/utilities/get_promoter_seqs.R")

# ------------------------------------------------------------------------------
# import data

# list of differentially expressed genes
deg <- read_csv(args$deg_csv) %>% 
  # use nominal p-value cuz AME will correct for p-values as well (don't want
  # to overcorrect)
  filter(
    p_val < 0.05 & if (args$logfc_pos) {
      avg_log2FC > 0
    } else {
      avg_log2FC < 0
    }
  )

# import GENCODE annotations to convert genes to genome ranges/coordinates
gencode <- rtracklayer::import(args$gencode_gtf) %>% 
  as.data.frame(.) %>% 
  dplyr::filter(type == "gene")

# set MEME suite database
meme_db <- args$meme_db

# ------------------------------------------------------------------------------
# get gene promoter sequences

# specify the genes
gencode <- dplyr::filter(gencode, gene_name %in% deg$gene)

message(sprintf("Getting promoter sequences for %i genes", nrow(gencode)))
promoter_seqs <- get_promoter_seqs(
  gencode_df = gencode,
  upstream = 1000,
  downstream = 100,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# ------------------------------------------------------------------------------
# run memes (AME, analysis of motif enrichment)

message("Running AME")
enr_motifs <- runAme(promoter_seqs, outdir = out_dir, database = meme_db)

message("Saving enriched motifs")
write_csv(x = enr_motifs, file = file.path(out_dir, "enriched_motifs.csv"))

# ------------------------------------------------------------------------------
# make heatmap

.plt <- plot_ame_heatmap(enr_motifs[1:30, ], value = -log10(adj.pvalue)) + 
  # need to increase margin cuz x-axis label is cut off
  theme(plot.margin = margin(t = 12, r = 12, b = 12, l = 24))

message("Saving heatmap")
ggsave(
  "ame_heatmap.png",
  plot = .plt,
  path = out_dir,
  width = 20,
  height = 5,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------

print(sessionInfo())
