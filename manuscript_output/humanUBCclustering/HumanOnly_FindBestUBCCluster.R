
rm(list=ls())

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/"
dt <- format(Sys.Date(),"%y%m%d")

inFile <- sprintf("%s/RLlineage_only/241107/UBC_withClusterAssignments_241107.qs",outDir)

seppMetadata <-  "/home/rstudio/isilon/private/projects/HumanMouseUBC/Sepp2024_metadata.txt"

outRoot <- sprintf("%s/UBCclusters_241112",outDir)
#outRoot <- sprintf("%s/UBCclusters_241118",outDir) #,dt)
if (!file.exists(outRoot)) dir.create(outRoot)

source("cluster_integrity_fromLeo.R")
dtm <- format(Sys.time(),"%y%m%d_%H%M")
logFile <- sprintf("%s/clusterIntegrity_%s.log",outRoot,dtm)

srat <- qs::qread(inFile)
con <- file(logFile)
#sink(con,split=TRUE)
#sink(con,type="message",append=TRUE)

tryCatch({
    # get age and sex for Sepp data
    cat("* Assigning age and sex columns...\n")
    mdata <- srat@meta.data
    age <- mdata$age
    sex <- mdata$sex
    sepp_age <- paste(sub(" wpc", "", mdata$Stage),"PCW")
    age[grep("Sepp",mdata$dataset_name)] <- sepp_age[grep("Sepp",mdata$dataset_name)]

    cat("\t* Assigning age to object...\n")
    srat$age <- age

    # extract sex information from Sepp samples
    seppFull <- read.delim(seppMetadata,sep="\t",h=T,as.is=T)
    seppSex2 <- seppFull[,c("orig.ident","TissueID")]
    seppSex2 <- seppSex2[!duplicated(seppSex2),]
    seppSex2 <- na.omit(seppSex2)
    spos <- gregexpr(" ",seppSex2$TissueID)
    seppSex2$sex <- sapply(1:length(spos),function(i) substr(seppSex2[i,2],
        spos[[i]][2]+1, spos[[i]][3]-1))
    
    cat("\t* MANUALLY ADDING SEX for SN297 (Male)\n")
    seppSex2$sex[which(seppSex2$orig.ident =="SN297")] <- "M"
    
    outFile <- sprintf("%s/SeppMetadata_Sex.txt",dirname(seppMetadata))
    write.table(seppSex2,file=outFile,sep="\t",col=T,row=F,quote=FALSE)

    idx <- grep("Sepp",mdata$dataset_name)
    matched_sex <- seppSex2$sex[
        match(mdata$orig.ident[idx],
        seppSex2$orig.ident)]
    mdata$sex[idx] <- matched_sex

    cat("\t* Assigning sex to object...\n")
    srat$sex <- mdata$sex

   # cat("\n* FILTERING FOR CELLS with non-zero EOMES expression\n")
   # eomes <- srat[["RNA"]]$counts["EOMES",]
   # srat$XPR_EOMES <- eomes > 0
   # srat <- subset(srat, XPR_EOMES == TRUE)
   # cat(sprintf("BEFORE: %i cells; AFTER: %i cells\n",
   #     length(eomes), dim(srat)[2]))

    cat(sprintf("\nSample breakdown:\n"))
    cat(sprintf("Num unique samples = %i\n", 
        length(unique(srat$orig.ident))))

    cat(sprintf("\nSex breakdown\n"))
    print(table(srat$sex,useNA="always"))
    cat(sprintf("\nAge breakdown\n"))
    print(summary(srat$age))
    cat(sprintf("\nDataset breakdown\n"))
    print(table(srat$dataset_name,useNA="always"))

    gencodeFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/gencode.v44.basic.annotation.gtf"

    message("reading gencode")
    gencode <- rtracklayer::import(gencodeFile)
    gencode <- as.data.frame(gencode)
    gencode <- subset(gencode,
        type %in% "gene" & gene_type %in% "protein_coding")

    profileClusters(
        dataset_srat=srat, 
        outDir=outRoot, 
        runSCT=FALSE, runDEG=FALSE,
        resSet=0.2, #c(0.4,0.6,0.8),
        dotplot_genes2plot=c("MKI67","WLS", "EOMES", "LMX1A", "OTX2", "PAX6","RBFOX3", "ATOH1","CALR","SOX4","SOX11","CRX","EYS","NRL"), 
        featureplot_genes2plot=c("EOMES","EYS","FOXP2","SOX4","CALR"),
        origClusterColumn="curtype",
        heatmap_topAnnot_column=c("age","sex"), 
        heatmap_log2FC_cutoffs=c(1,1.5,2),
        heatmap_addGenes=c("SOX4","SOX11","EYS","CRX","NRL"),
        gencodeList=gencode
)
}, error=function(ex){
    print(ex)
}, finally={
    cat("Closing log files.")
 #   sink(NULL)
 #   sink(NULL)
})

