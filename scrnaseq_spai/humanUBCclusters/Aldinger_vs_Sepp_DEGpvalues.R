# plot nominal p for UBC1 genes in Aldinger vs Sepp
rm(list=ls())
require(ggplot2)
require(ggrepel)


#DEGresults <- "/home/rstudio/isilon/private/llau/results/integrated/20240715/all_tested_genes.csv"

DEGresults <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBC1_marker_analysis/Integrated_RL_Hsonly/DEG_notfiltered_241029.txt"

outDir <- "/home/rstudio/isilon/private/projects/HumanMouseUBC/UBC1_marker_analysis"

dt <- format(Sys.Date(),"%y%m%d")

orig_dat <- read.delim(DEGresults,sep="\t",h=T,as.is=T)

for  (cl in unique(orig_dat$cluster)) {
message(cl)
dat <- subset(orig_dat,cluster %in% cl)
dat <- subset(dat, 
    Aldinger_full_cerebellum_human_p_val < 0.05 & 
    Sepp_full_cerebellum_human_p_val < 0.05
)

dat <- dat[,c(1,2,6,7,14)]
dat$avgLogFC <- (dat[,2]+dat[,4])/2
dat$Ald <- dat$Aldinger_full_cerebellum_human_avg_log2FC
dat$Sepp <- dat$Sepp_full_cerebellum_human_avg_log2FC

#dat2 <- subset(dat, abs(avgLogFC) >1.5)

p <- ggplot(dat, aes(x=Ald,y=Sepp,label=feature)) + 
    geom_point(aes(colour=avgLogFC)) +
    scale_colour_gradient2()
#p <- p + geom_label_repel(
#   data=dat2, 
#    aes(x=Ald, y=Sepp, label=feature),
#    max.overlaps=80, 
#    size=2) #+ xlim(c(0,1.5)) + ylim(c(0,2))

cor <- cor.test(dat$Ald,dat$Sepp)

p <- p + ggtitle(sprintf("%s DEG - Sepp vs Aldinger (Pearson = %1.2f, p<%1.2e)",   cl, cor$estimate, cor$p.value))
p <- p + xlab(sprintf("log2FC(%s vs other), Aldinger 2021",cl))
p <- p + ylab(sprintf("log2FC(%s vs other), Sepp 2023",cl))
p <- p + geom_hline(yintercept=0,linetype="dotted")
p <- p + geom_vline(xintercept=0,linetype="dotted")
#p <- p + theme(text=element_text(size=20))

ggsave(p, 
    file=sprintf("%s/AldvsSepp_%s_DEGnotfiltered_%s.pdf",outDir,cl,dt)
)

}