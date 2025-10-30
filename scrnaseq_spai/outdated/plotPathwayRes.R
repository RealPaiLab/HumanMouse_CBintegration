rm(list=ls())

inDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/results/integrated_HsFullCB/20241031/cca/UBCclusters_241112/0.2_res"

pFiles <-  dir(inDir, pattern="pathwayEnrichment.csv")


pathList <- list()
for (f in pFiles){
    dat <- read.delim(sprintf("%s/%s",inDir,f),sep=",",h=T,as.is=T)
    dat <- subset(dat, p_value < 0.05)
    cat(sprintf("%s: %i pathways Q < 0.05\n", 
        f,nrow(dat)))
    dat$ptrans <- -log10(dat$p_value)
    print(dat$term_name)
    #p <- ggplot(dat,aes(y=ptrans)) + geom_bar()
    #p <- p + geom_col()
    #ggsave(p,file=sprintf("%s/pathwayBarplot.png",inDir)
    cat("\n\n")
    pathList[[f]] <- dat$term_name

}