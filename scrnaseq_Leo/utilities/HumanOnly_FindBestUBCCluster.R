rm(list=ls())

outDir <- "/.mounts/labs/pailab/private/llau/results/human/20241125"
dt <- format(Sys.Date(),"%y%m%d")

inFile <- "/.mounts/labs/pailab/private/llau/results/human/20241108/rl/0.4_ald_ubc.qs"

outRoot <- sprintf("%s/UBCclusters_241125",outDir) #,dt)
if (!file.exists(outRoot)) dir.create(outRoot)

source("/u/llau/software/mb_scrnaseq/MB_scRNAseq/scrnaseq_Leo/utilities/cluster_integrity_fromLeo.R")
dtm <- format(Sys.time(),"%y%m%d_%H%M")
logFile <- sprintf("%s/clusterIntegrity_%s.log",outRoot,dtm)

srat <- qs::qread(inFile)
#con <- file(logFile)
sink(logFile,split=TRUE)
#sink(logFile,type="message"==,append=TRUE)

tryCatch({
    # get age and sex for Sepp data
    cat("* Assigning age and sex columns...\n")
    mdata <- srat@meta.data
    srat$dataset_name <- "Aldinger"
    age <- mdata$age
    sex <- mdata$sex

    cat("\n* FILTERING FOR CELLS with non-zero EOMES expression\n")
    eomes <- srat[["RNA"]]@counts["EOMES",]
    srat$XPR_EOMES <- eomes > 0
    srat <- subset(srat, XPR_EOMES == TRUE)
    cat(sprintf("BEFORE: %i cells; AFTER: %i cells\n",
        length(eomes), dim(srat)[2]))

    cat(sprintf("\nSample breakdown:\n"))
    cat(sprintf("Num unique samples = %i\n", 
        length(unique(srat$orig.ident))))

    cat(sprintf("\nSex breakdown\n"))
    print(table(srat$sex,useNA="always"))
    cat(sprintf("\nAge breakdown\n"))
    print(summary(srat$age))
    cat(sprintf("\nDataset breakdown\n"))
    print(table(srat$dataset_name,useNA="always"))

    #gencodeFile <- "/home/rstudio/isilon/private/projects/FetalHindbrain/anno/gencode.v44.basic.annotation.gtf"

    #message("reading gencode")
    #gencode <- rtracklayer::import(gencodeFile)
    #gencode <- as.data.frame(gencode)
    #gencode <- subset(gencode,
        #type %in% "gene" & gene_type %in% "protein_coding")

    profileClusters(
        dataset_srat=srat, 
        outDir=outRoot, 
        runSCT=TRUE, runDEG=TRUE,
        resSet=c(0.2,0.4,0.6,0.8),
        dotplot_genes2plot=c("MKI67","WLS", "EOMES", "LMX1A", "OTX2", "PAX6","RBFOX3", "ATOH1","CALR","SOX4","SOX11","CRX","EYS","NRL"), 
        featureplot_genes2plot=c("EOMES","EYS","FOXP2","SOX4","CALR"),
        origClusterColumn="snn_res.0.4",
        heatmap_topAnnot_column=c("age","sex"), 
        heatmap_log2FC_cutoffs=c(1,1.5,2),
        heatmap_addGenes=c("SOX4","SOX11","EYS","CRX","NRL")
        #gencodeList=gencode
)
}, error=function(ex){
    print(ex)
}, finally={
    cat("Closing log files.")
    sink()
 #   sink(NULL)
})