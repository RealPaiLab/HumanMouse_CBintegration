### UTILS ###
#' get Cavalli metadata 
get_cavalli_meta <- function(cavalli_meta_file) {
  cavalli_meta <- readxl::read_xlsx(cavalli_meta_file, 
                                    sheet = 1, skip = 1, n_max = 765,
                                    col_names = T, trim_ws = T, na = "NA")
  cavalli_meta <- cavalli_meta[,1:10]
  cavalli_meta$Subgroup <- factor(cavalli_meta$Subgroup, 
                                  levels = c("WNT", "SHH", "Group3", "Group4")
  )
  
  return(cavalli_meta)
}

#' get Cavalli expression matrix
get_cavalli_expr <- function(cavalli_expr_file) {
  cavalli_expr <- read.csv(cavalli_expr_file, header = T, sep = "\t")
  print(dim(cavalli_expr))
  
  # transpose to sample X gene
  expr_tmp <- t(cavalli_expr)
  colnames(expr_tmp) <- expr_tmp["HGNC_symbol_from_ensemblv77",]
  expr_tmp <- expr_tmp[6:nrow(expr_tmp),]
  expr_tmp <- as.data.frame(expr_tmp)
  expr_tmp$Study_ID <- rownames(expr_tmp)
  
  return(expr_tmp)
}

#' get Hendrikse RNA-seq DESeq2 dds
get_hendrikse_dds <- function(hendrikse_meta_file, hendrikse_expr_file) {
  require(DESeq2)
  # meta
  hendrikse_meta <- readxl::read_xlsx(hendrikse_meta_file, 
                                      sheet = "Group_3_4_Bulk_Samples"
  )
  col_data <- as.data.frame(hendrikse_meta)
  col_data$Subgroup <- as.factor(col_data$Subgroup)
  rownames(col_data) <- col_data$Sample_ID
  col_data$Sample_ID <- NULL
  
  # expr
  mtx <- as.matrix(read.table(hendrikse_expr_file, check.names = F))
  
  # two samples were not matched to meta data, typo? 
  # MDT-AP-2109_MB-2109_Primary_A18561 -> MDT-AP-0721_MB-2109? 
  # MDT-AP-2120_MB-2120_Primary_A25296 -> MDT-AP-2532_MB-2120?
  # MDT-AP-1145_MB-1145_Primary_A08576 ~ the only shh tumour
  # remove for accuracy
  mtx <- mtx[,! colnames(mtx) %in% c("MDT-AP-2109_MB-2109_Primary_A18561", 
                                     "MDT-AP-2120_MB-2120_Primary_A25296", 
                                     "MDT-AP-1145_MB-1145_Primary_A08576")] 
  
  common_name <-substr(colnames(mtx), 1, 19)
  col_data <- col_data[common_name,]
  colnames(mtx) <- common_name
  
  dds <- DESeqDataSetFromMatrix(countData = mtx,
                                colData = col_data,
                                design= ~ Subgroup)
  dds <- DESeq(dds, parallel = T)
  
  return(dds)
}

#' plot a red X with message
#' @param e (character) The message
plot_error <- function(e) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.1, label = e, color = "black", size = 2) +
    geom_segment(aes(x = 0.2, y = 0.8, xend = 0.8, yend = 0.2), color = "red", size = 2) +  # Diagonal line
    geom_segment(aes(x = 0.2, y = 0.2, xend = 0.8, yend = 0.8), color = "red", size = 2) +  # Other diagonal line
    ggtitle("Error occurred!") +
    theme_void()
  
  return(p)
}

#' plot patched figures of mut, cavalli, and hendrikse RNA
#' @param gene (character) Name of the interested gene
#' @param mutation (GRanges) Mutations with a meta col named as "class"
#' @param gencode (GRanges) GENCODE full info
#' @param commonSeqLevels (logical) Default TRUE to keep only common seqLevels
#' @param meta_df (data.frame) Meta data df
#' @param expr (data.frmae) Cavalli expression count 
#' @param comparison (list) A list of paired subgroup for wilcox test. e.g. list(c(1,2))
#' @param dds (DESeqDataSet) DDS object of DESeq2 on hendrikse count matrix
#' @param subgroups (character) A vector of subgroups to keep for survival plot
plot_gene_sum <- function(gene, 
                          mutation, gencode, commonSeqLevels = T,
                          meta_df, expr, comparison = NULL,
                          dds, 
                          subgroups = c("Group3", "Group4")
                          ) {
  require(patchwork)
  
  message(sprintf("--- Plotting summary for gene: %s", gene))
  
  tryCatch(
    p1 <- plot_mutation(gene, mutation = mutation, gencode = gencode),
    error = function(e) {
      p1 <<- plot_error(e)
    }
  )
  
  tryCatch(
    p2 <- plot_bulkRNA_cavalli(gene, meta_df = meta_df, expr = expr),
    error = function(e) {
      p2 <<- plot_error(e)
    }
  )
  
  tryCatch(
    p3 <- plot_bulkRNA_hendrikse(gene, dds = dds),
    error = function(e) {
      p3 <<- plot_error(e)
    }
  )
  
  tryCatch({
    l <- plot_survival(gene, meta_df = meta_df, expr = expr, subgroups = subgroups)
    count <- sum(l$data.survtable[l$data.survtable$time == 0, "n.risk"])
    
    surv_p <- l$plot + 
      labs(title = sprintf("%s (Array; n = %i)", paste(subgroups, collapse = " & "), count)) + 
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
    table_p <- l$table + 
      theme(plot.title = element_blank(),legend.position = "none")

    p4 <- (surv_p/table_p) + plot_layout(heights = c(3,1))
  },
    error = function(e) {
      p4 <<- plot_error(e)
    }
  )
  
  p <- ((p2/p3)|p1|p4) + plot_layout(widths = c(0.38, 0.22, 0.38))
  
  return(p)
}
