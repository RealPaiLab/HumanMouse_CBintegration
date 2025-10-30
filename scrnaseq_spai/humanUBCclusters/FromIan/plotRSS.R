# read PySCENIC rss output and plot top regulons in Aldinger RL SVZ cluster
rm(list=ls())
library(ggplot2)
library(ggrepel)

inFile <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/20250319/human_ubcs.rss.csv"
outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/integrated_human_ubc/PySCENIC"

dt <- format(Sys.time(), "%y%m%d")
outDir <- sprintf("%s/%s", outDir, dt)
if (!dir.exists(outDir)) {
    dir.create(outDir, recursive=FALSE)
}

logFile <- sprintf("%s/log.txt", outDir)

sink(logFile,split=TRUE)
tryCatch({
dat <- read.delim(inFile,sep=",",h=T,as.is=T)

useOrder <- paste("UBC_",c(7,3,8,2,6,4,1,0,5),sep="")

for (cl in useOrder) {
    cat(cl)
    cat("\n")
    cur <- subset(dat, X == cl)[,-1]
    nm <- names(cur)
    cur <- as.numeric(cur)
    idx <- order(cur,decreasing=TRUE)
    cur <- cur[idx]; nm <- nm[idx]
    
    cur <- as.data.frame(cur)
    cur$nm <- nm
    colnames(cur) <- c("score","regulon")

    print(head(cur,8))
    cat("\n")
    print(cur[1:8,2,drop=F])
    cat("\n\n")

    maxi <- min(20,nrow(cur))
    p <- ggplot(cur[1:maxi,], 
        aes(x=reorder(regulon, score, decreasing=TRUE),y=score)) + 
        geom_point(size=2) + 
        geom_text_repel(
            data=cur[1:8,], 
            aes(label=regulon), col="red",
            max.overlaps=10,
            size=8,
            hjust=1) +
        ylab("Regulon Specificity Score") + 
        xlab("Decreasing rank order") + 
        theme_minimal() + 
        theme(  text=element_text(size=30),
                axis.text.y=element_text(size=24),
                axis.text.x = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black")
        )
        ggtitle(sprintf("Aldinger PySCENIC: Cluster %s", cl))
    ggsave(p, file=sprintf("%s/cluster_%s.pdf", outDir, gsub(" ",".",cl)), 
        width=6, height=7, dpi=300,unit="in")

}

},error=function(e){
    cat("Error: ", e$message, "\n")
}, finally={
    cat("Done.\n")
    sink()
})

