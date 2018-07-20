#setwd("/Users/lpsmith/R/ASCAT_test")
args <- commandArgs(trailingOnly = TRUE)
p = "999"
if (length(args) > 0) {
  p=args[1]
}
gamma <- 240
if(length(args) > 1){
  gamma <- as.integer(args[2])
}
print ("Using a gamma of:")
print(gamma)
.libPaths(c("~/R/", .libPaths()))
#library("devtools")
#install("~/R/ls_copynumber")
library("lscopynumber")
lucianTest()
sapply(list.files(pattern="[.]r$", path="~/R/ls_copynumber/", full.names=TRUE), source);

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
logR_T <- logR_T[which(lna < (ncol(logR_T) - 2)),]
logR_T$Chr <- as.numeric(logR_T$Chr)
logR_T$Position <- as.numeric(logR_T$Position)
logR_T <- logR_T[order(logR_T$Chr, partial=logR_T$Position),]
im <- as.data.frame(apply(logR_T, 2, as.numeric))
im <- winsorize(im, gamma=40)
#im <- imputeMissing(im, "pcf")
save(im, file=sprintf("%s_im.RData", p))
logsegs <- multipcf(im, gamma=gamma, w=samplew)
save(logsegs, file=sprintf("%s_logsegs.RData", p))

BAF_T[BAF_T[,1] == "X", 1] <- 23
BAF_T[BAF_T[,1] == "Y", 1] <- 24
for (i in colnames(BAF_N)[-c(1:2)]) {
  BAF_T[which(BAF_N[,i] > 0.75 | BAF_N[,i] < 0.25),i] <- NA
}
nna <- apply(BAF_T[,-c(1,2)], 1, function(x){length(which(is.na(x)))})
baf <- BAF_T[which(nna < length(BAF_N)-2),]
baf$Chr <-  as.numeric(baf$Chr)
baf[,-c(1,2)] <- 0.5+abs(baf[,-c(1,2)]-0.5)
baf <- winsorize(baf, gamma=40)
#baf <- imputeMissing(baf, "pcf")
bafsegs <- multipcf(baf, gamma=gamma, w=samplew)
save(bafsegs, file=sprintf("%s_bafsegs.RData", p))
bafsegs <- bafsegs[which(bafsegs$n.probes > 10),]

blsegs <- rbind(logsegs, bafsegs)
chrom.all <- c()
start.all <- c()
end.all <- c()
for(i in 1:22){
  s.chrom <- subset(blsegs, as.numeric(blsegs[, "chrom"]) == i)
  start <- unique(as.numeric(s.chrom[, "start.pos"]))
  end <-  unique(as.numeric(s.chrom[, "end.pos"]))
  end1 <- 0
  start.final <- c()
  end.final <- c()
  count <- 1
  while (end1 < max(start) & max(end) > min(start[start > end1])) {
    start1 <- min(start[start > end1])
    end1 <- min(end[end > start1])
    start.final <- c(start.final, start1)
    end.final <- c(end.final, end1)
    count <- count + 1    
  }
  chrom.all <- c(chrom.all, rep(i, length(start.final)))
  start.all <- c(start.all, start.final)
  end.all <- c(end.all, end.final)
}

mat.pos <- cbind(as.numeric(chrom.all), as.numeric(start.all), as.numeric(end.all))
colnames(mat.pos) <- c("Chr", "Start", "End")

nlog <- nbaf <- c()
for (chr in unique(mat.pos[,"Chr"])) {
  subPos <- mat.pos[which(mat.pos[,"Chr"] == chr),,drop=F]
  subLog <- logR_T[which(logR_T[,"Chr"] == chr),,drop=F]
  subBaf <- baf[which(baf[,"Chr"] == chr),,drop=F]
  nlog <- c(nlog, apply(subPos, 1, function(coord,pos){length(which(pos >= coord["Start"] & pos <= coord["End"]))},
                        pos=subLog[,"Position"]))
  nbaf <- c(nbaf, apply(subPos, 1, function(coord,pos){length(which(pos >= coord["Start"] & pos <= coord["End"]))},
                        pos=subBaf[,"Position"]))
}
mat.pos <- mat.pos[which(nbaf>0),]
nlog <- nlog[which(nbaf>0)]
nbaf <- nbaf[which(nbaf>0)]
cnSegs <- cbind(mat.pos, nlog, nbaf)
colnames(cnSegs) <- c(colnames(mat.pos), "nLogR", "nBAF")
write.table(cnSegs, file=sprintf("%s_copynumber_segments.txt", p), sep="\t", row.names=F, quote=F)

source("aspcf.R")
getLogBaf <- function(j, logR, Baf) {
  ilogR <- logR[,j]
  ibaf <- Baf[,j]
  yi2 <- ibaf[which(!is.na(ibaf))]
  ##sd1 <- getMad(ilogR)
  allBflip <- yi2
  allBflip[yi2 > 0.5] <- 1 - yi2[yi2 > 0.5]

  sd2 <- getMad(allBflip)
  if (is.na(sd2))
    sd2 <- 0
  ## Center data around zero (by subtracting 0.5) and estimate mean
  if(length(yi2)== 0){
    mu <- 0
  }else{
    mu <- mean(abs(yi2-0.5))
  }
  
  ## Make a (slightly arbitrary) decision concerning branches
  ## This may be improved by a test of equal variances
  if(sqrt(sd2^2+mu^2) < 2*sd2){
    mu <- 0
  }
  return(c(mean(ilogR,na.rm=T), mu+0.5))
}

logBafPerSeg <- function(pos, subLog, subBaf) {
  subLog2 <- subLog[which(subLog[,1] >= pos[1] & subLog[,1] <= pos[2]),]
  subBaf2 <- subBaf[which(subBaf[,1] >= pos[1] & subBaf[,1] <= pos[2]),]
  return(lapply(2:ncol(subLog2), getLogBaf, subLog2, subBaf2))
}

segLog <- segBaf <- segCN <- matrix(NA,nc=(ncol(logR_T)-2), nr=nrow(mat.pos))
colnames(segLog) <- colnames(segBaf) <- colnames(logR_T)[-c(1:2)]
for (chr in unique(mat.pos[,"Chr"])) {
  print(chr)
  subPos <- mat.pos[which(mat.pos[,"Chr"] == chr),-1,drop=F]
  subLog <- logR_T[which(logR_T[,"Chr"] == chr),-1,drop=F]
  subBaf <- BAF_T[which(BAF_T[,"Chr"] == chr),-1,drop=F]
  
  logBaf <- matrix(unlist(apply(subPos, 1, logBafPerSeg, subLog, subBaf)), nc=((ncol(logR_T)-2)*2), byrow=T)
  segLog[which(mat.pos[,"Chr"] == chr),] <- logBaf[,seq(1, ncol(logBaf)-1, 2)]
  segBaf[which(mat.pos[,"Chr"] == chr),] <- logBaf[,seq(2, ncol(logBaf), 2)]
}
segLog <- cbind(mat.pos, segLog)
segBaf <- cbind(mat.pos, segBaf)
save(segLog, file=sprintf("%s_noSmear_segLog.RData",p))
save(segBaf, file=sprintf("%s_noSmear_segBaf.RData",p))
