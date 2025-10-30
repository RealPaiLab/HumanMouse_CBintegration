rm(list=ls())

inFile <- "/home/rstudio/isilon/private/nghafouri/OriginsOfMB/data/Erickson_tumour/250805/erickson_metadata.txt"



dat <- read.delim(inFile, header=T, sep=" ")

selID <- c("HSMB0008", "HSMB0058", "HSMB0061",
    "HSMB0006", "HSMB0074","HSMB0226","HSMB0048",
    "HSMB0054","HSMB0181",
    "HSMB0135","HSMB0002","HSMB0130")

 tmp <- subset(dat, sample_id %in% selID)
blah <- tmp[!duplicated(tmp$sample_id),]
print(table(blah$subtype))