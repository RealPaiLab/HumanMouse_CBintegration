# ==============================================================================
# Given a gene (e.g., BRCA1), count number of mutations that are present in MB
# stratified by location (e.g. exon, promoter...) and tumour subtype.
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
    "--gencode_gtf", "/isilon/CBL_scRNAseq-archived/data/src/gencode/GRCh37/gencode.v19.annotation.gtf.gz",
    "--out_dir", "/CBL_scRNAseq/results/tumour/20240225"
  ))
} else {
  arg_list <- parser$parse_args()
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
    "/isilon/CBL_scRNAseq-archived/data/src/ICGC/PBCA-DE/simple_somatic_mutation.open.PBCA-DE.tsv.gz",
    "/isilon/CBL_scRNAseq-archived/data/src/ICGC/PEME-CA/simple_somatic_mutation.open.PEME-CA.tsv.gz"
  ),
  .f = \(x) {
    message(sprintf("\n***Reading in %s***", x))
    m <- read_tsv_chunked(x, DataFrameCallback$new(subset_mb)) %>% 
      dplyr::select(c(
        icgc_donor_id,
        project_code,
        submitted_sample_id,
        chromosome:mutation_type,
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

# >>> import patient metadata
donor_meta <- readxl::read_xlsx(
  "/isilon/CBL_scRNAseq-archived/data/src/medulloblastoma-genomics/genetic-variants/Northcott_2017/41586_2017_BFnature22973_MOESM1_ESM.xlsx",
  sheet = 1
)
# <<<

# ------------------------------------------------------------------------------
# identify coordinates of a gene

filter_range <- function(gencode_gr, gene_name, seq_type) {
  return(gencode_gr[gencode_gr$gene_name == gene_name & gencode_gr$type == seq_type, ])
}

get_gene_ranges <- function(gencode, gene_name) {
  # get range of gene body (transcription start to end)
  gene_body <- filter_range(
    gencode_gr = gencode,
    gene_name = gene_name,
    seq_type = "gene"
  ) %>% 
    reduce() %>% 
    granges()
  
  # get ranges of gene exons
  gene_exon <- filter_range(
    gencode_gr = gencode,
    gene_name = gene_name,
    seq_type = "exon"
  ) %>% 
    reduce() %>% 
    granges()
  
  # get range of gene promoter region
  gene_promoter <- promoters(gene_body, upstream = 1000, downstream = 0)
  
  # make GRangesList
  grlist <- GRangesList(
    "gene_body" = gene_body,
    "exon" = gene_exon,
    "promoter" = gene_promoter
  )
  
  return(grlist)
}

brca1 <- get_gene_ranges(gencode, "BRCA1")

# ------------------------------------------------------------------------------
# find mutations in the gene regions (gene body, exon, promoter)

get_mutation_gene_overlaps <- function(
  gene_region_ranges,
  mb_mutation_ranges
) {
  overlaps <- map(
    .x = gene_region_ranges,
    .f = \(gr) {
      olp <- GenomicRanges::findOverlaps(
        query = mb_mutation_ranges,
        subject = gr,
        minoverlap = 1L,
        type = "any",
        select = "all",
        ignore.strand = TRUE
      )
      
      mut_df <- mb_mutation_ranges[from(olp), ] %>% 
        as.data.frame()
      
      return(mut_df)
    }
  ) %>% 
    list_rbind(names_to = "seq_type")
  
  return(overlaps)
}

# ------------------------------------------------------------------------------
# join tumour subtype (from patient metadata) to the mutation

add_mb_subgroup <- function(
  mutation_gene_overlaps,
  donor_meta,
  by = c("submitted_sample_id" = "PID")
) {
  overlaps <- left_join(
    x = mutation_gene_overlaps,
    y = donor_meta,
    by = by
  ) %>%
    dplyr::rename(subgroup = SUBGROUP) %>%
    mutate(subgroup = factor(.$subgroup, levels = c("WNT", "SHH", "Group 3", "Group 4")))
  
  return(overlaps)
}

# ------------------------------------------------------------------------------
# get number of mutations for a list of genes

genes <- c("BRCA1", "CBFA2T2", "OTX2", "MYC", "MYCN", "CTNNB1", "CNTNAP2")

walk(
  .x = genes,
  .f = \(gene) {
    # get gene coordinates
    gr <- get_gene_ranges(gencode = gencode, gene_name = gene)
    
    # find mutations in gene
    overlaps <- get_mutation_gene_overlaps(
      gene_region_ranges = gr,
      mb_mutation_ranges = mb_mutations
    )
    
    # add tumour subtype data
    overlaps <- add_mb_subgroup(
      mutation_gene_overlaps = overlaps,
      donor_meta = donor_meta
    )
    
    write_csv(
      x = overlaps,
      file = file.path(arg_list$out_dir, paste0(gene, ".csv"))
    )
    
    plt <- ggplot(overlaps, aes(x = seq_type, fill = subgroup)) + 
      geom_bar(position = position_dodge()) + 
      labs(title = gene, x = "gene region", y = "number of mutations") + 
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
      theme_classic() + 
      theme(axis.ticks = element_line(colour = "black"),
            axis.text = element_text(colour = "black"))
    
    ggsave(
      filename = paste0(gene, ".png"),
      plot = plt,
      path = arg_list$out_dir,
      width = 8,
      height = 4,
      units = "in",
      dpi = 600
    )
  }
)

# ------------------------------------------------------------------------------

message("\n***SESSION INFO***\n")
print(sessionInfo())

