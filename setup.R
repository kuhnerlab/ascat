library("devtools")
setwd("/Users/lpsmith/R/ASCAT_test")
install("copynumber")
library("copynumber")
lucianTest()

#args <- commandArgs(trailingOnly = TRUE)
#p=args[1]
p = "702"
sapply(list.files(pattern="[.]r$", path="copynumber/", full.names=TRUE), source);

## reformat data
## fill in paths to the data files here!

workingdir <- getwd()

logR_T <- read.table(sprintf(paste(workingdir, "/%s_logR.txt", sep=""), p),check.names=FALSE,stringsAsFactors=FALSE)
BAF_T <- read.table(sprintf(paste(workingdir, "/%s_BAF.txt", sep=""), p),check.names=FALSE,stringsAsFactors=FALSE)
BAF_N <- read.table(sprintf(paste(workingdir, "/%s_Normal_BAF.txt", sep=""), p),check.names=FALSE,stringsAsFactors=FALSE)

samplew <- rep(1, ncol(logR_T)-2)
samplew[grep("whole", colnames(logR_T)[-c(1:2)])] <- 2

library(limma)
logR_T[logR_T[,1] == "X",1] <- 23
logR_T[logR_T[,1] == "Y",1] <- 24
lna <- apply(logR_T[,-c(1,2)], 1, function(x){length(which(is.na(x)))})
logR_T <- logR_T[which(lna < (ncol(logR_T) - 2)/2),]
logR_T$Chr <- as.numeric(logR_T$Chr)
logR_T$Position <- as.numeric(logR_T$Position)
logR_T <- logR_T[order(logR_T$Chr, partial=logR_T$Position),]
im <- as.data.frame(apply(logR_T, 2, as.numeric))
im <- winsorize(im, gamma=40)
#im <- imputeMissing(im, "pcf")
save(im, file=sprintf("%s_im.RData", p))
save(samplew, file=sprintf("%s_samplew.RData", p))


#install("copynumber")
#logsegs <- multipcf(im, gamma=40, w=samplew)
