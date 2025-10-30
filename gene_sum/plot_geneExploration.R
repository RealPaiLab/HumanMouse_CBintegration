rm(list=ls())

library(ggplot2)
library(dplyr)
library(ggpubr)


### FUNCTIONS ###
## UTILS ##
# TODO: borrowed from previous scripts with brutal modifications, can improve
#' Load GENCODE gene Transcription Start Sites or gene body
#' @param genocode (GRanges) GENOCDE info
#' @param gene_type (character) Either "all" or "protein_coding" genes
#' @param region (character) Either "tss" or "gene_body"
#' @return (GRanges) Return a Granges object of the gene TSS
get_gencode_anno <- function(gencode,
                             gene_type = "all", 
                             region = "tss") {
  ### Filter gene type ###
  genes <- as.data.frame(gencode)
  if (tolower(gene_type) == "all") {
    genes <- subset(genes, type == "gene")
  } else if (tolower(gene_type) == "protein_coding") {
    genes <- subset(genes, gene_type %in% "protein_coding" & type == "gene")
  } else {
    stop(sprintf("Wrong gene_type provided: %s", gene_type))
  }
  
  ### Set ranges ###
  if (tolower(region) == "tss") {
    genes$TSS <- genes$start
    # tss flip for reverse strand
    genes$TSS[which(genes$strand=="-")] <- genes$end[which(genes$strand=="-")] 
    genes$l <- genes$TSS
    genes$r <- genes$TSS
  } else if (tolower(region) == "gene_body") {
    genes$l <- genes$start
    genes$r <- genes$end
  } else {
    stop(sprintf("Wrong region provided: %s", region))
  }
  
  geneGR <- GRanges(
    genes$seqnames,
    IRanges(genes$l, genes$r),
    name=genes$gene_name, 
    strand = genes$strand
  ) 
  
  ### Output ###
  return(geneGR)
}


#' get upstream tss region
#' @param tss (GRanges) TSS
#' @param up_len (numeric) Number of base pairs up of TSSs. Default 1kb
#' @return (GRanges)
get_promoters <- function(tss, up_len = 1000) {
  ### set promoter gr ###
  promoter <- tss
  loc <- start(promoter)
  strand <- as.character(strand(tss))
  if (strand == "+") {
    start(promoter) <- loc - up_len
    end(promoter) <- loc - 1
  } else if (strand == "-") {
    start(promoter) <- loc + 1
    end(promoter) <- loc + up_len
  } else {
    stop("Missing or wrong TSS strand info")
  }
  
  return(promoter)
}


#' find muts in a gene
#' @param gene (character) Name of the interested gene
#' @param mutation (GRanges) Mutations with a meta col named as "class"
#' @param gencode (GRanges) GENCODE full info
#' @param commonSeqLevels (logical) Default TRUE to keep only common seqLevels
#' @return (GRanges) overlap mutations in a gene; if not in gencode, return NULL
find_ol_muts <- function(gene, mutation, gencode, 
                         commonSeqLevels = T) {
  # QC #
  if (! gene %in% gencode$gene_name) {
    print(sprintf("Gene %s is not included in the provided GENCODE.", gene))
    return(NULL)
  }
  
  if (commonSeqLevels) {
    commonSeqLevels <- paste0("chr", c(1:22, "X", "Y", "M"))
    mutation <- mutation[seqnames(mutation) %in% commonSeqLevels]
    seqlevels(mutation) <- commonSeqLevels 
  }
  
  # prep gencode intervals and types #
  sub_gencode <- gencode[gencode$gene_name == gene & 
                           gencode$type %in% c("gene", "exon")]
  uptss <- get_promoters(tss = get_gencode_anno(gencode = sub_gencode))
  sub_gencode <- sub_gencode[order(sub_gencode$type != "exon")]
  sub_gencode <- sub_gencode[! duplicated(sub_gencode)]
  
  # append upstream TSS region to the gencode
  uptss$name <- NULL
  uptss$type <- "upTSS"
  sub_gencode <- c(sub_gencode, uptss)
  
  # aggregate muts #
  olp <- findOverlapPairs(mutation, sub_gencode)
  ol_mut <- olp@first
  ol_mut$type <- olp@second$type
  
  return(ol_mut)
}


#' Generate mutation count matrix of genes in the given gene regions
#' @param genes (character) vector of interested genes
#' @param groups (character) vector of interested MB groups (Grp3, Grp4, SHH, WNT)
#' @param types (character) vector of gene region types (upTSS, exon, gene)
#' @return (matrix) matrix with genes in row, group_type in col
generate_mut_df <- function(genes, 
                            groups =  c("Grp3", "Grp4"),
                            types = c("upTSS","exon", "gene"), ...) {
  # prepare mtx to fill
  mtx <- matrix(nrow = length(genes), ncol = length(groups)*length(types))
  rownames(mtx) <- genes
  colnames(mtx) <- expand.grid(groups, types) %>% 
    arrange(factor(Var1, levels = c("Grp4", "Grp3", "SHH", "WNT")), 
            factor(Var2, levels = c("upTSS","exon", "gene"))) %>%
    mutate(res = paste(Var1, Var2, sep = "_")) %>%
    pull(res)
  
  for (gene in genes) {
    ol_mut <- find_ol_muts(gene = gene, 
                           mutation = combined_mut, 
                           gencode = gencode)
    if (is.null(ol_mut)) {
      mtx[gene,] <- NA
      next
    }
    
    for (group in groups) {
      sum_df <- ol_mut %>% 
        as.data.frame() %>%
        arrange(factor(type, levels = c("upTSS","exon", "gene"))) %>%
        filter(class == group) %>%
        group_by(type) %>%
        summarize(N_mut = length(unique(name)))
      
      for (type in types) {
        v <- ifelse(type %in% sum_df$type, 
                    sum_df[sum_df$type == type,]$N_mut,
                    0
        )
        mtx[gene, paste(group, type, sep = "_")] <- v
      }
    }
  }
  
  # format colnames
  colnames(mtx) <- sub("_gene", "_gene-other", colnames(mtx))
  colnames(mtx) <- sub("_upTSS", "_upstream-1kb", colnames(mtx))
  
  return(mtx)
}


## PLOT ##
#' plot mutation profile of a gene by classification
#' @param gene (character) Name of the interested gene
#' @param mutation (GRanges) Mutations with a meta col named as "class"
#' @param gencode (GRanges) GENCODE full info
#' @param commonSeqLevels (logical) Default TRUE to keep only common seqLevels
#' @return (list) a list contains plot and summary data frame
plot_mutation <- function(gene, mutation, gencode, 
                          commonSeqLevels = T) {
  ol_mut <- find_ol_muts(gene, mutation, gencode, commonSeqLevels)
  
  sum_df <- ol_mut %>% 
    as.data.frame() %>%
    arrange(factor(type, levels = c("upTSS","exon", "gene"))) %>%
    distinct(class, name, .keep_all = TRUE) %>% 
    group_by(class, type) %>%
    summarize(N_mut = length(unique(name)))
    
  tumour_count <- table(unique(mcols(mutation))$class)
  
  # plot #
  p <- ggplot(data = sum_df, aes(x = class, y = N_mut, fill = type)) + 
    geom_bar(stat = "identity", position = "stack") + 
    scale_fill_manual(values = c(gene = "#4D4D4D", 
                                 exon = "#56B4E9",
                                 upTSS = "#009E73"
                                 ), 
                      labels = c(gene = "Gene (exons excluded)", 
                                 exon = "Exon",
                                 upTSS = "1kb Upstream TSS"
                                 )) +
    labs(y = "Number of mutated tumours", 
         title = sprintf("%s (PCAWG; n = %i)", gene, sum(tumour_count)), 
         fill = NULL, 
         subtitle = paste(names(tumour_count), tumour_count, 
                          sep = ": ", collapse = "; ")
         ) +
    geom_text(aes(label = N_mut, group = type), 
              position = position_stack(vjust = 0.5),
              size = 5, color = "white") +
    theme_minimal() + 
    theme(text = element_text(size = 14), 
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_blank(),
          legend.position = "top")
    
  return(list(plot = p, sum_df = sum_df))
}


#' plot bulk RNA from cavalli dataset
#' @param gene (character) The target gene
#' @param meta_df (data.frame) Meta data df
#' @param expr (data.frmae) Cavalli expression count 
#' @param comparison (list) A list of paired subgroup for wilcox test. e.g. list(c(1,2))
#' @return (gg)
plot_bulkRNA_cavalli <- function(gene, meta_df, expr, comparison = NULL) {
  # set comparison pairs for test
  if (is.null(comparison)) {
    comparison <- list(c("Group3", "Group4"),
                       c("Group4", "SHH"),
                       c("SHH", "WNT"),
                       c("Group3", "SHH"), 
                       c("Group3", "WNT"), 
                       c("Group4", "WNT")
                       )
  }
  
  # QC #
  problematic_genes <- c("CCDC7", "KLK9", "WDR92") # gene names are duplicated in expr
  
  if (! gene %in% colnames(expr)) {
    stop(sprintf("Gene %s is not included in Cavalli data set.", gene))
  } else if ( gene %in% problematic_genes) {
    stop(sprintf("Gene %s is not allowed.", gene))
  }
  
  # merge expr and meta #
  merged <- expr %>% 
    select("Study_ID", gene) %>% 
    mutate_at(.vars = gene, .funs = as.numeric) %>%
    inner_join(meta_df, by = "Study_ID")
  
  # plot #
  p <- merged %>% ggplot(aes_string(x = "Subgroup", y = gene, color = "Subgroup")) + 
    geom_violin(show.legend = F) + 
    geom_jitter(alpha = 0.2, show.legend = F) +
    scale_color_manual(values = c(WNT = "#557aab", 
                                  SHH = "#c65960",
                                  Group3 = "#ebdb94",
                                  Group4 = "#629083"
                                  )) +
    labs(y = "normalized log2 expr", 
         title = sprintf("%s (Array; n = %i)", gene, nrow(merged))) + 
    stat_compare_means(
      method = "wilcox.test",
      comparisons = comparison,
      label = "p.format",
      label.x = 1.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -40),
          text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.x = element_blank(),
          legend.position = "top"
          ) + 
    stat_summary(fun.data = function(x) {
      return(data.frame(y = min(x) - 0.4, label = paste("n=", length(x), sep = "")))
    }, geom = "label", color = "black")
  
  return(p)
}


#' plot bulk RNA from hendrikse dataset
#' @param gene (character) The target gene
#' @param dds (DESeqDataSet) DDS object of DESeq2 on hendrikse count matrix
#' @return (gg)
plot_bulkRNA_hendrikse <- function(gene, dds) {
  # package installation guide: 
  # BiocManager::install(c("AnnotationFilter", "ensembldb", "ComplexHeatmap"))
  # install.packages("RNAseqQC")
  require(RNAseqQC)
  dds_gene_symbols <- stringr::str_split(rownames(dds), "___", simplify = T)[,1]
  
  if (! gene %in% dds_gene_symbols) {
    stop(sprintf("Gene %s is not included in Hendrikse data set.", gene))
  }
  
  id <- rownames(dds)[dds_gene_symbols == gene]
  
  p <- plot_gene(id, dds, x_var = "Subgroup", show_plot = F)$data %>% 
    ggplot(aes(x = Subgroup, y = count, color = Subgroup)) + 
    geom_violin(show.legend = F) + 
    geom_jitter(alpha = 0.2, show.legend = F) +
    scale_color_manual(values = c(Group3 = "#ebdb94",Group4 = "#629083")) +
    labs(y = "log2(1 + normalized count)", 
         title = sprintf("%s (RNA-seq; n = %i)", gene, ncol(dds))) + 
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("Group3", "Group4")),
      label = "p.format",
      label.x = 1.5,
    ) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -40),
          text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.x = element_blank(),
          legend.position = "top"
    ) + 
    stat_summary(fun.data = function(x) {
      return(data.frame(y = min(x) - 0.4, label = paste("n=", length(x), sep = "")))
    }, geom = "label", color = "black")
  
  return(p)
}


#' plot survival plot of a gene
#' @param gene (character) Name of the interested gene
#' @param meta_df (data.frame) Meta data df
#' @param expr (data.frmae) Cavalli expression count 
#' @param subgroups (character) A vector of subgroups to keep
#' @return (list)
plot_survival <- function(gene, meta_df, expr,
                          subgroups = c("Group3", "Group4")) {
  require(survival)
  require(survminer)
  
  # QC #
  problematic_genes <- c("CCDC7", "KLK9", "WDR92") # gene names are duplicated in expr
  
  if (! gene %in% colnames(expr)) {
    stop(sprintf("Gene %s is not included in Cavalli data set.", gene))
  } else if ( gene %in% problematic_genes) {
    stop(sprintf("Gene %s is not allowed.", gene))
  }
  
  # meta
  sub <- meta_df %>%
    filter(Subgroup %in% subgroups) %>%
    select(Study_ID, Age, Gender, 
           Dead, `Met status (1 Met, 0 M0)`, `OS (years)`,
           Subgroup) %>%
    filter(!is.na(Dead)) %>%
    mutate_at(.vars = c("Age", "Dead", "Met status (1 Met, 0 M0)", "OS (years)"), 
              .funs = as.numeric) %>%
    mutate(`OS (months)` = `OS (years)`*12)
  
  # categorize by gene expr and combine with meta
  expr %>% 
    select(Study_ID, all_of(gene)) %>% 
    mutate_at(.vars = gene, .funs = as.numeric) %>%
    mutate(median = median(!!sym(gene))) %>% 
    mutate(actual_median = median(as.numeric(expr$SOX4))) %>%
    head()
  merged <- expr %>% 
    select(Study_ID, all_of(gene)) %>% 
    mutate_at(.vars = gene, .funs = as.numeric) %>%
    inner_join(sub, by = "Study_ID") %>%
    mutate(category = ifelse(!!sym(gene) > median(!!sym(gene), na.rm = T), 
                             sprintf("%s High", gene),
                             sprintf("%s Low", gene)
    ))
  
  # survival curve
  fit <- survfit(Surv(`OS (months)`, Dead) ~ category, data = merged)
  l <- ggsurvplot(fit, data = merged, 
                  conf.int = TRUE, pval = TRUE, 
                  risk.table = TRUE,
                  risk.table.col = "strata",
                  legend.title = "",
                  xlab = "Time (months)",
                  palette = c("#E69F00", "#56B4E9"),
                  legend.labs = c(sprintf("%s High", gene), sprintf("%s Low", gene)),
                  pval.method = T
  )
  
  return(l)
}



# plot g34 mut pyramid
#' plot Group 3/4 MB mutation pyramimd (left snv_indel, right snv_indel_sv)
#' @param cluster (charater) The target cluster in dict. e.g. "UBC1"
#' @param dict (list) genes of cluster. e.g. list(UBC1 = ubc1_up_tf, UBC2 = ubc2_up_tf)
#' @param combined_mut_snvindel, 
#' @param combined_mut,
#' @param gencode
#' @return (ggplot)
plot_mut_pyramid <- function(cluster, dict,
                             combined_mut_snvindel, 
                             combined_mut_snvindelsv,
                             gencode
                             ) {
  # g4
  grp4_df <- NULL
  for (tf in dict[[cluster]]) {
    grp <- c("Grp4")
    # snv-indel
    tmp <- plot_mutation(tf, 
                         combined_mut_snvindel[combined_mut_snvindel$class %in% grp], 
                         gencode)$sum_df
    tmp <- tmp[tmp$class %in% grp,]
    tmp$TF <- tf
    tmp$mutType <- "snv-indel"
    grp4_df <- rbind(grp4_df, tmp)
    
    # snv-indel-sv
    tmp <- plot_mutation(tf, 
                         combined_mut[combined_mut$class %in% grp], 
                         gencode)$sum_df
    tmp <- tmp[tmp$class %in% grp,]
    tmp$TF <- tf
    tmp$mutType <- "snv-indel-sv"
    grp4_df <- rbind(grp4_df, tmp)
  }
  
  miss <- dict[[cluster]][!dict[[cluster]] %in% grp4_df$TF]
  if (length(miss) != 0) {
    grp4_df <- rbind(grp4_df, 
                     data.frame(class = "Grp4", 
                                type = "gene", 
                                N_mut = 0, 
                                TF = miss, 
                                mutType = "snv-indel-sv")
    )
  }
  
  grp4_tf_order <- grp4_df %>%
    filter(mutType == "snv-indel-sv") %>%
    group_by(TF) %>%
    summarise(sum = sum(N_mut), .groups = "drop") %>%
    arrange(desc(sum)) %>%
    pull(TF)
  
  grp4_muts_plot <- ggplot(data = grp4_df, aes(x = factor(TF, levels = grp4_tf_order), 
                                               y = ifelse(mutType == "snv-indel", -N_mut, N_mut), 
                                               fill = type)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c(gene = "#4D4D4D", 
                                 exon = "#56B4E9",
                                 upTSS = "#009E73"
    ), 
    labels = c(gene = "Gene (exons excluded)", 
               exon = "Exon",
               upTSS = "1kb Upstream TSS"
    )) +
    scale_y_continuous(labels = abs) +
    coord_flip() + 
    labs(fill = NULL) +
    ylab("Number of mutated tumours") +
    xlab(sprintf("%s Upregulated TFs", cluster)) + 
    ggtitle(sprintf("%s Upregulated TFs Mutation Counts (%s)", 
                    cluster,
                    paste(grp, collapse = "&")
    ), 
    subtitle = sprintf(
      "1kb Upstream TSS > Exon > Gene (exons excluded)\nn=%i",
      length(unique(combined_mut[combined_mut$class %in% grp]$name)))
    ) +
    annotate("text", x = length(grp4_tf_order)+1, y = 2, label = "SNV-INDEL-SV", 
             hjust = 0, vjust = 1, size = 5, fontface = "bold", color = "grey50") +  
    annotate("text", x = length(grp4_tf_order)+1, y = -2, label = "SNV-INDEL", 
             hjust = 1, vjust = 1, size = 5, fontface = "bold", color = "grey50") +
    geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 1) +
    geom_text(aes(label = N_mut, group = type), 
              position = position_stack(vjust = 0.5),
              size = 5, color = "white") +
    theme_minimal() + 
    theme(text = element_text(size = 14), 
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.position = "top")
  
  
  # g3
  grp3_df <- NULL
  for (tf in dict[[cluster]]) {
    grp <- c("Grp3")
    # snv-indel
    tmp <- plot_mutation(tf, 
                         combined_mut_snvindel[combined_mut_snvindel$class %in% grp], 
                         gencode)$sum_df
    tmp <- tmp[tmp$class %in% grp,]
    tmp$TF <- tf
    tmp$mutType <- "snv-indel"
    grp3_df <- rbind(grp3_df, tmp)
    
    # snv-indel-sv
    tmp <- plot_mutation(tf, 
                         combined_mut[combined_mut$class %in% grp], 
                         gencode)$sum_df
    tmp <- tmp[tmp$class %in% grp,]
    tmp$TF <- tf
    tmp$mutType <- "snv-indel-sv"
    grp3_df <- rbind(grp3_df, tmp)
  }
  miss <- dict[[cluster]][!dict[[cluster]] %in% grp3_df$TF]
  if (length(miss) != 0) {
    grp3_df <- rbind(grp3_df, 
                     data.frame(class = "Grp3", 
                                type = "gene", 
                                N_mut = 0, 
                                TF = miss, 
                                mutType = "snv-indel-sv")
    )
  }
  
  grp3_tf_order <- grp3_df %>%
    filter(mutType == "snv-indel-sv") %>%
    group_by(TF) %>%
    summarise(sum = sum(N_mut), .groups = "drop") %>%
    arrange(desc(sum)) %>%
    pull(TF)
  
  grp3_muts_plot <- ggplot(data = grp3_df, aes(x = factor(TF, levels = grp4_tf_order), 
                                               y = ifelse(mutType == "snv-indel", -N_mut, N_mut), 
                                               fill = type)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c(gene = "#4D4D4D", 
                                 exon = "#56B4E9",
                                 upTSS = "#009E73"
    ), 
    labels = c(gene = "Gene (exons excluded)", 
               exon = "Exon",
               upTSS = "1kb Upstream TSS"
    )) +
    scale_y_continuous(labels = abs) +
    coord_flip() + 
    labs(fill = NULL) +
    ylab("Number of mutated tumours") +
    xlab(sprintf("%s Upregulated TFs", cluster)) + 
    ggtitle(sprintf("%s Upregulated TFs Mutation Counts (%s)",
                    cluster,
                    paste(grp, collapse = "&")
    ), 
    subtitle = sprintf(
      "1kb Upstream TSS > Exon > Gene (exons excluded)\nn=%i",
      length(unique(combined_mut[combined_mut$class %in% grp]$name)))
    ) +
    annotate("text", x = length(grp3_tf_order)+1, y = 2, label = "SNV-INDEL-SV", 
             hjust = 0, vjust = 1, size = 5, fontface = "bold", color = "grey50") +  
    annotate("text", x = length(grp3_tf_order)+1, y = -2, label = "SNV-INDEL", 
             hjust = 1, vjust = 1, size = 5, fontface = "bold", color = "grey50") +
    geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 1) +
    geom_text(aes(label = N_mut, group = type), 
              position = position_stack(vjust = 0.5),
              size = 5, color = "white") +
    theme_minimal() + 
    theme(text = element_text(size = 14), 
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.position = "top")
  
  # sum plot
  require(patchwork)
  return(grp4_muts_plot | grp3_muts_plot)
}









  

  
  
  
