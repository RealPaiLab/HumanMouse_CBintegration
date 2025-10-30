# ==============================================================================
# Take ICGC mutations and output a file that can be used as input to the Ensembl
# Variant Effect Predictor (VEP).
# ==============================================================================

library(argparse)
library(tidyverse)
library(GenomicRanges)

# parse command line arguments
parser <- ArgumentParser()

parser$add_argument(
  # GENCODE GTF annotation file
  "--gencode_gtf",
  default = NULL,
  required = TRUE
)
parser$add_argument(
  # output directory
  "--out_dir",
  default = NULL,
  required = TRUE
)

if (interactive()) {
  arg_list <- parser$parse_args(c(
    "--gencode_gtf", "/.mounts/labs/pailab/src/gencode/GRCh37/gencode.v19.annotation.gtf.gz",
    "--out_dir", "/.mounts/labs/pailab/private/icheong/CBL_scRNAseq/results/tumour/ICGC/20250226"
  ))
} else {
  arg_list <- parser$parse_args()
}

message(sprintf("Saving files to %s", getwd()))
if (!dir.exists(arg_list$out_dir)) {
  dir.create(arg_list$out_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# load data

# gencode file (genome annotation)
gencode <- rtracklayer::import(arg_list$gencode_gtf)

# >>> mutation data

# filter for MB samples only
subset_mb <- function(x, pos) {
  x <- dplyr::filter(x, grepl("MDT-AP|ICGC_MB", submitted_sample_id))
  return(x)
}

mb_mutations <- map(
  .x = c(
    "/.mounts/labs/pailab/src/ICGC/PBCA-DE/simple_somatic_mutation.open.PBCA-DE.tsv.gz",
    "/.mounts/labs/pailab/src/ICGC/PEME-CA/simple_somatic_mutation.open.PEME-CA.tsv.gz"
  ),
  .f = \(x) {
    message(sprintf("\n***Reading in %s***", x))
    m <- read_tsv_chunked(x, DataFrameCallback$new(subset_mb)) %>% 
      dplyr::select(c(
        icgc_donor_id,
        project_code,
        submitted_sample_id,
        chromosome:mutated_to_allele,
        gene_affected
      )) %>% 
      mutate(
        chromosome = paste0("chr", chromosome)
      ) %>% 
      # remove duplicate rows (same mutation may affect multiple transcripts)
      unique
    return(m)
  }
) %>% 
  list_rbind() %>% 
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    start.field = "chromosome_start",
    end.field = "chromosome_end"
  )
# <<<

# ------------------------------------------------------------------------------
# identify coordinates of a gene

filter_range <- function(gencode_gr, gene_name, seq_type) {
  return(gencode_gr[gencode_gr$gene_name == gene_name & gencode_gr$type == seq_type, ])
}

get_gene_body <- function(gencode, gene_name) {
  # get range of gene body (transcription start to end)
  gr <- filter_range(
    gencode_gr = gencode,
    gene_name = gene_name,
    seq_type = "gene"
  ) %>% 
    reduce() %>% 
    granges()

  return(gr)
}

# ------------------------------------------------------------------------------
# find mutations in the gene body

get_mutation_gene_overlaps <- function(
  gene_region_ranges,
  mb_mutation_ranges
) {
  overlaps <- GenomicRanges::findOverlaps(
    query = mb_mutation_ranges,
    subject = gene_region_ranges,
    minoverlap = 1L,
    type = "any",
    select = "all",
    ignore.strand = TRUE
  )

  mut_df <- mb_mutation_ranges[from(overlaps), ] %>% 
    as.data.frame()

  return(mut_df)
}

# ------------------------------------------------------------------------------
# save coordinates for BRCA1 mutations in default VEP input format

# get gene body coordinates
gr <- get_gene_body(gencode = gencode, gene_name = "BRCA1")

# find mutations
overlaps <- get_mutation_gene_overlaps(
  gene_region_ranges = gr,
  mb_mutation_ranges = mb_mutations
)

# save output for VEP
# (see https://www.ensembl.org/info/docs/tools/vep/vep_formats.html)
overlaps %>%
  mutate(
    seqnames = str_remove(string = seqnames, pattern = "chr"),
    start = case_when(
      str_detect(string = mutation_type, pattern = "insertion") ~ end + 1,
      .default = start
    ),
    allele = sprintf("%s/%s", reference_genome_allele, mutated_to_allele),
    strand = "+"
  ) %>%
  select(seqnames, start, end, allele, strand) %>%
  # add a positive control (from Ander Diaz-Navarro; PMID 26200345 fig 2)
  add_row(
    seqnames = "9",
    start = 139390152,
    end = 139390152,
    allele = "T/C",
    strand = "+"
  ) %>%
  # save to file
  write.table(
    file = file.path(arg_list$out_dir, "BRCA1_vep_input.txt"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

