### Code from Shraddha's GitHub

# Prioritize CRISPRi candidates
# rm(list=ls())

library(readxl)
library(gplots) # venn intersection items
# library(ggplot2)
# library(ggrepel)

# FetalHindbrain_Epigenetics/anno/CRISPR_prioritization
inDir <- "./criteria_for_filtering"

geneSets <- list(
  vz_svz_TF="Aldinger_TFs/vz_svz/ame.tsv",
  svz_eubc_TF="Aldinger_TFs/svz_eubc/ame.tsv",
  RL_EGL_DNAm="RL_EGL_DNAm/RL_EGL_DMR_ProteinCodingGenes_ame_fromAlexChan.tsv",
  EpiDB="EpiDB/Proteins.csv",
  MB_alterations="MB_alterations/MB_Genes.xlsx",
  MB_cancerindex="MB_alterations/MBgenes_cancerindexorg_230814.txt",
  D425_CRISPR="D425_CRISPR/SD2ori.xls",
  depmap_xpr="depmap/D425/Expression_Public_22Q4_subsetted.csv",
  depmap_crispr="depmap/D425/CRISPR_Depmap_Public_22Q4_plus_Score_Chronos.txt",
  humanTF="HumanTFs_CCBR/TF_names_v_1.01.txt",
  UBCTF="Fetalbrain_snRNAseq/HumanSpecificUBC_Genes_TFs.txt",
  CRISPR_essentials="D425_CRISPR/CRISPR_common_essentials.csv"
)

# dt <- format(Sys.Date(),"%y%m%d")
# outDir <- sprintf("%s/prioritize_%s", inDir, dt)
# if (!file.exists(outDir)) dir.create(outDir)
outDir <- "./criteria_for_filtering_results"

geneList <- list()

# AME results
for (nm in c("vz_svz_TF","svz_eubc_TF")){
  message(nm)
  dat <- read.delim(sprintf("%s/%s",inDir,geneSets[[nm]]))
  x <- substr(dat$motif_ID,1,regexpr("_",dat$motif_ID)-1) 
  geneList[[nm]] <- cbind(gene=unique(x),type=nm)
  message(sprintf("\t%i genes", length(unique(x))))
}

# keep TFs that are upregulated in the stem cells
message("Filter TFs enriched in SVZ-eUBC & upregulated in SVZ")
dat <- read.delim(sprintf("%s/Aldinger_TFs/svz_eubc.txt",inDir),
                  sep="\t",h=T,as.is=T)
upreg <- subset(dat, p_val_adj<0.05 & avg_log2FC > 0)
message(sprintf("\t%i upregulated genes", nrow(upreg)))
geneList[["svz_eubc_TF"]] <- intersect(
  geneList[["svz_eubc_TF"]][,1],rownames(upreg))
geneList[["svz_eubc_TF"]] <- union(geneList[["svz_eubc_TF"]], rownames(upreg))

message(sprintf("\t%i SVZ-UBC TFs upregulated in SVZ + SVZ upreg genes", 
                length(geneList[["svz_eubc_TF"]])))
message(sprintf("{ %s }", 
                paste(geneList[["svz_eubc_TF"]],collapse=","))
)
rm(dat,upreg)

message("Filter TFs enriched in VZ-SVZ & upregulated in SVZ")
dat <- read.delim(sprintf("%s/Aldinger_TFs/vz_svz.txt",inDir),
                  sep="\t",h=T,as.is=T)
upreg <- subset(dat, p_val_adj<0.05 & avg_log2FC > 0)
message(sprintf("\t%i upregulated genes", nrow(upreg)))
geneList[["vz_svz_TF"]] <- intersect(
  geneList[["vz_svz_TF"]][,1],rownames(upreg))
geneList[["vz_svz_TF"]] <- union(geneList[["vz_svz_TF"]], rownames(upreg))
message(sprintf("\t%i VZ-SVZ TFs also upregulated in VZ", 
                length(geneList[["vz_svz_TF"]])))
##message(sprintf("{ %s }", 
##    paste(geneList[["vz_svz_TF"]],collapse=","))
##)

dt <- format(Sys.Date(),"%y%m%d")
write.table(geneList[["vz_svz_TF"]], 
            file=sprintf("%s/vz_svz_Upreg_vz_%s.txt",outDir, dt))
write.table(geneList[["svz_eubc_TF"]], 
            file=sprintf("%s/svz_eubc_Upreg_svz_%s.txt",outDir, dt))

epi <- read.delim(sprintf("%s/%s",inDir,geneSets$EpiDB),sep=",")
geneList[["EpiDB"]] <- unique(epi$HGNC_symbol)

mb <- read_excel(sprintf("%s/%s",inDir,geneSets$MB_alterations),skip=3)
x <- as.data.frame(mb)
mb2 <- read.delim(sprintf("%s/%s",inDir,geneSets$MB_cancerindex),sep="\t",h=T,as.is=T)
geneList[["MB_MutOE"]] <- unique(x[,1],mb2$Gene)

genesOnly <- geneList

xs <- venn(genesOnly, show.plot=FALSE)
print( lengths(attributes(xs)$intersections))

message("vz-svz TF + EpiDB + MB MutOE")
print(sort(attributes(xs)$intersections[["vz_svz_TF:EpiDB:MB_MutOE"]]))

message("svz-eUBC TF + MB MutOE")
g1 <- attributes(xs)$intersections[["svz_eubc_TF:MB_MutOE"]]
print(sort(g1))

message("vz-svz TF + MB MutOE")
g2 <- attributes(xs)$intersections[["vz_svz_TF:MB_MutOE"]]
print(sort(g2))

message("svz-eUBC TF + EpiDB")
g3 <- attributes(xs)$intersections[["svz_eubc_TF:EpiDB"]]
print(sort(g3))

message("svz-eUBC TF + MB MutOE")
g4 <- attributes(xs)$intersections[["svz_eubc_TF:MB_MutOE"]]
print(sort(g4))

xpr <- read.delim(sprintf("%s/%s",inDir,geneSets$depmap_xpr),sep=",")
xpr <- t(xpr); 
xpr <- data.frame(
  gene=rownames(xpr)[-(1:6)], 
  xpr=as.numeric(xpr[-(1:6)])
)
crispr <- read.delim(sprintf("%s/%s", inDir,geneSets$depmap_crispr),sep=",")
crispr <- t(crispr)
crispr <- data.frame(
  gene=rownames(crispr)[-(1:6)],
  crispr=as.numeric(crispr[-(1:6)])
)

tf <- read.delim(sprintf("%s/%s", inDir,geneSets$humanTF),header=FALSE)[,1]
ubctf <- read.delim(sprintf("%s/%s",inDir,geneSets$UBCTF), header=FALSE)[,1]

selectedGenes <- unique(c(geneList[["vz_svz_TF"]], geneList[["svz_eubc_TF"]]))
message(sprintf("From scRNAseq analysis = %i targets",length(selectedGenes)))

selepi <- intersect(selectedGenes, geneList$EpiDB)
selTF <- intersect(selectedGenes, tf)
selectedGenes <- union(selepi,selTF)
message(sprintf("--> (Is Epi) or (IsTF) = %i targets",length(selectedGenes)))

ess <- read.delim(sprintf("%s/%s", inDir, geneSets$CRISPR_essentials),
                  sep=" ")
ess <- rownames(ess)    
selectedGenes <- setdiff(selectedGenes, ess)
message(sprintf("Exclude CRISPR common essential genes -> %i targets",
                length(selectedGenes)))

# validate on Depmap
# https://depmap.org/portal/interactive/?filter=&regressionLine=false&associationTable=false&x=slice%2FChronos_Combined%2F7794%2Fentity_id&y=slice%2Fexpression%2F7794%2Fentity_id&color= 
d425 <- merge(x=xpr,y=crispr,by="gene")
d425 <- subset(d425, gene %in% selectedGenes)
message(sprintf("-->In D425 dataset: %i targets", nrow(d425)))


# set yes/no flags for additional gene attributes
d425$MB_mutOE <- d425$gene %in% geneList$MB_MutOE
d425$svz_eUBC_TF <- d425$gene %in% geneList$svz_eubc_TF
d425$IsEpiReg <- d425$gene %in% geneList$EpiDB
d425$IsTF <- d425$gene %in% tf
d425$IsUBCTF <- d425$gene %in% ubctf

d425 <- d425[order(d425$crispr, -d425$xpr),]

d425$selLabel <- rep("", nrow(d425))
idx <- union(1:10,
             #which(d425$xpr > 1 & d425$crispr < -0.25),
             #which(d425$gene %in% c("FOXP2","FOXP1","ZIC1")),
             which(d425$MB_mutOE == TRUE))
#which(d425$gene %in% ubctf))
d425$selLabel[idx] <- d425$gene[idx]
#d425 <- subset(d425, IsTF==TRUE | IsEpiReg==TRUE) 



p <- ggplot(d425,aes(x=crispr,y=xpr)) + 
  geom_point(aes(col=!IsEpiReg)) + 
  geom_vline(xintercept=0,col="grey75",linetype=2) +
  geom_hline(yintercept=d425$xpr[50],col="grey75",linetype=2) + 
  geom_text_repel(
    aes(label = selLabel),
    min.segment.length = 0.1
  ) + 
  labs(
    title = "D425 cells (G3 MB cell line)",
    x = "CRISPR gene dependency score",
    y = "log2(TPM) expression"
  ) + 
  scale_colour_manual(
    values = unname(pals::trubetskoy()[4:5]),
    name = "Is epigenetic\nregulator?",
    labels = c("Yes", "No")
  ) + 
  # scale_colour_brewer(
  #   palette = "Dark2",
  #   # rev = 
  #   name = "Is epigenetic\nregulator?",
  #   labels = c("Yes", "No")
  # ) + 
  theme_classic() + 
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = c(0.85, 0.85), 
    legend.box.background = element_rect(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
ggsave(
  "crispr_depmap.png",
  plot = p,
  width = 4,
  height = 5,
  units = "in",
  dpi = 1200
)
# outFile <- sprintf("%s/candidateGenes_D425view_%s.png",outDir,dt)
# ggsave(p2,file=outFile)
#print(p2)
#dev.off()

###geom_text(label=d425$selLabel, nudge_y=0.2,size=5, )
