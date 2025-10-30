

inFile <- "/home/rstudio/isilon/private/projects/MB_multiome/input/AkdesHarmanci/finalChrMat_morestringent_v2.rda"
load(inFile)

# count cells with malignancies
x <- colSums(abs(finalChrMat_morestringent))
normal <- names(which(x < 1))
write.table(normal, file="normal_cells.txt",col=F,row=F,quote=F)
