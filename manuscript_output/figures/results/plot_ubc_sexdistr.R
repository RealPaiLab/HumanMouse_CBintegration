# plot sex distribution for UBC 1 and 2
rm(list=ls())
library(Seurat)
library(ggplot2)
library(dplyr)

source("/home/rstudio/isilon/private/icheong/CBL_scRNAseq/software/utilities/cluster_barplot.R")

ubcFile <- "/home/rstudio/isilon/private/icheong/CBL_scRNAseq/results/integrated/20240825/ubc_subset.qs"
seppMetadataFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/Sepp2024_metadata.txt"

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/sexDistr"
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/%s", outDir,dt)

# ---------------------------------------

#' Custom `ggplot2` theme based on `theme_classic()`.
#'
#' @param base_size Base font size in points
#' @param base_family Base font family
#' @param base_line_size Base size for line elements
#' @param base_rect_size Bsae size for rect elements
#'
#' @return A `ggplot2` theme.
#'
theme_classic2 <- function(
  base_size = 11,
  base_family = "",
  base_line_size = base_size / 22,
  base_rect_size = base_size / 22
) {
  theme_classic(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace% 
    theme(
      axis.text = element_text(colour = "black", size = rel(0.8)),
      axis.ticks = element_line(colour = "black"),
      legend.key.size = unit(0.8, "lines"),
      legend.box.spacing = unit(0.25 * base_size, "points"),
      plot.title = element_text(size = rel(1.2), hjust = 0.5),
      plot.tag = element_text(face = "bold", hjust = 0.5, vjust = 0.5),
      plot.tag.location = "plot",
      strip.background = element_blank(),
      strip.text = element_text(size = rel(1)),
      complete = TRUE
    )
}


# -------------------
if (!dir.exists(outDir)){ dir.create(outDir, recursive=FALSE)}
logFile <- sprintf("%s/log.txt",outDir)
sink(logFile, split=TRUE)
tryCatch({
cat("Reading human-mouse UBC\n")
srat <- qs::qread(ubcFile)

# subset for species == human
srat <- subset(srat, subset=species == "human")
mdata <- srat[[]]
cat(sprintf("Have %i human UBC cells\n\n",nrow(mdata) ))

sex <- mdata$sex
cat("Adding Sepp sample sex info\n")
seppSexInfo <- sprintf("%s/SeppMetadata_Sex.txt",dirname(seppMetadataFile))
seppSex <- read.table(seppSexInfo, header=TRUE, sep="\t", stringsAsFactors=FALSE)

sex <- mdata[,c("TissueID","sex","dataset_name")]
midx <- match(sex$TissueID, seppSex$TissueID)
if (all.equal(seppSex$TissueID[midx],sex$TissueID)!=TRUE){
    cat("tissueID didn't match\n")
}
idx <- grep("Sepp",sex$dataset_name)
sex[idx,2] <- seppSex$sex[midx[idx]]

if (all.equal(rownames(mdata), rownames(sex))!=TRUE){
    cat("mdata and sex rownames didn't match\n")
}
srat$combined_sex <- sex$sex

cat(sprintf("After merging, sex distribution\n"))
md <- srat[[]]
print(table(md[,c("dataset_name","combined_sex")]),
useNA="always")

cat("\nBy proportion\n")
# print proportion of sex per dataset
for (d in unique(md$dataset_name)){
    tab <- table(md$combined_sex[md$dataset_name==d])
    prop <- round(100*tab/sum(tab),2)
    cat(sprintf("Dataset: %s\n",d))
    print(prop, useNA="always")
}   

cat("Print counts of UBC subclusters by sex\n")
print(table(md[,c("subclusters","combined_sex")]), useNA="always")

cat("Both datasets combined: % male, by subcluster\n")
for (k in sort(unique(md$subclusters))){
    cur <- subset(md, subclusters == k)
    curm <- sum(cur$combined_sex %in% "M")
    
    cat(sprintf("%s: %% male = %1.2f%%\n",
        k, 100*(curm/nrow(cur))))
}

cat("\nJUST ALDINGER: % male, by subcluster\n")
md2 <- subset(md, dataset_name %in% "Aldinger_full_cerebellum_human")
for (k in sort(unique(md2$subclusters))){
    cur <- subset(md2, subclusters == k)
    curm <- sum(cur$combined_sex %in% "M")
    
    cat(sprintf("%s: %% male = %1.2f%%\n",
        k, 100*(curm/nrow(cur))))
}
cat("\n")

srat$clean_dataset_name <- gsub("_full_cerebellum_human","",md$dataset_name)
srat$subclusters <- sub("UBC_","UBC",srat$subclusters)

# for all cluster barplots, convert the y-axis labels to %
.plt1 <- cluster_barplot(
    srat,
    split.by = "combined_sex",
    group.by = "subclusters",
    position = "fill"
  ) + ggtitle("Sex by UBC subcluster - all datasets") + scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ theme_classic2() + 
    ylab("Proportion of cells")

suppressMessages(ggsave(.plt1, file=sprintf("%s/sex_by_subcluster_allUBC.pdf",outDir)))

.plt2 <- cluster_barplot(
    subset(srat, subset=dataset_name %in% c("Aldinger_full_cerebellum_human")),
    split.by = "combined_sex",
    group.by = "subclusters",
    position = "fill"
  ) +  ggtitle("Aldinger dataset") + scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ theme_classic2()

suppressMessages(ggsave(.plt2, file=sprintf("%s/sex_by_subcluster_AldingerOnly.pdf",outDir)))

.plt3 <- cluster_barplot(
    subset(srat, subset=dataset_name %in% c("Sepp_full_cerebellum_human")),
    split.by = "combined_sex",
    group.by = "subclusters",
    position = "fill"
  ) +  ggtitle("Sepp dataset") + scale_y_continuous(labels = scales::percent_format(accuracy = 1))+  theme_classic2()

suppressMessages(ggsave(.plt3, file=sprintf("%s/sex_by_subcluster_SeppOnly.pdf",outDir)))

.plt4 <- cluster_barplot(
    srat,
    split.by = "combined_sex",
    group.by = "clean_dataset_name",
    position = "fill"
  ) +  ggtitle("Sex by dataset") + scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ theme_classic2()
suppressMessages(ggsave(.plt4, file=sprintf("%s/sex_by_dataset.pdf",outDir)))

# write a table of counts by dataset, subcluster and sex. aggregate by subcluster and dataset
md <- srat[[]]
tab <- as.data.frame(table(md[,c("clean_dataset_name","subclusters","combined_sex")]))
tab <- tab[order(tab$subclusters, tab$clean_dataset_name, tab$combined_sex),]
colnames(tab) <- c("Dataset","UBC Subcluster","Sex","Cell Count")
outFile <- sprintf("%s/UBC_sex_by_dataset_subcluster_counts.txt",outDir)
write.table(tab, file=outFile, sep="\t", quote=FALSE, row.names=FALSE)

},error=function(ex){
    print(ex)
},finally={
    sink(NULL)
}
)


