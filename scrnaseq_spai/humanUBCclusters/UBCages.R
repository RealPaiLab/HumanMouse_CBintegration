library(Seurat)

UBC_Seurat <- "/home/rstudio/isilon/public/HumanMouseUBC/data/UBC.Harmony.RDS"

    cat("Reading UBC clusters ...")
    cat(sprintf("Input file:\n%s\n", UBC_Seurat))
    t0 <- Sys.time()
    srat <- readRDS(UBC_Seurat);
    cat("done in ", Sys.time()-t0, "\n")

    ages <- as.integer(sub(" PCW","",srat$age))
    srat$agewk <- ages
    mdata <- srat[[]]

