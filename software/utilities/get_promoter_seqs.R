#' Get promoter sequences for genes in a GENCODE dataframe format.
#'
#' @param gencode_df A dataframe containing GTF annotations. (Import from a GTF
#'   file and convert to a dataframe.)
#' @param upstream Number of bases upstream of the TSS
#' @param downstream Number of bases downstream of the TSS
#' @param genome BSgenome.Hsapiens.UCSC.hg38
#'
#' @return DNAStingSet object of promoter sequences for the given GENCODE
#'   dataframe
#' 
get_promoter_seqs <- function(
  gencode_df,
  upstream = 1000,
  downstream = 100,
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
) {
  # get the GRanges of the genes
  genes <- GRanges(
    seqnames = gencode_df$seqnames,
    ranges = IRanges(start = gencode_df$start,
                     end = gencode_df$end,
                     names = gencode_df$gene_name), 
    strand = gencode_df$strand,
    gene_name = gencode_df$gene_name
  )
  
  promoter_ranges <- IRanges::promoters(genes, upstream = upstream, downstream = downstream)
  
  promoter_seqs <- memes::get_sequence(
    regions = promoter_ranges,
    genome = genome
  )
  
  return(promoter_seqs)
}
