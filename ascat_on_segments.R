# load libraries
source("pierre_ascat.ko.mk.R")
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
p=args[1]
print(args[2])
doubled = args[2]

# we substitute reading raw data, hoping that will help
workingdir <- getwd()
logR_T <- read.table(sprintf(paste(workingdir, "/%s_logR.txt", sep=""),p),stringsAsFactors=F,check.names=F)
BAF_T <- read.table(sprintf(paste(workingdir, "/%s_BAF.txt", sep=""),p),stringsAsFactors=F,check.names=F)
# yes, we really read the same file twice here!  
logR_N <- read.table(sprintf(paste(workingdir, "/%s_logR.txt", sep=""),p),stringsAsFactors=F,check.names=F)
BAF_N <- read.table(sprintf(paste(workingdir, "/%s_Normal_BAF.txt", sep=""),p),stringsAsFactors=F,check.names=F)

# read file with start and end position of each segment, plus # SNPs
cnSegs <- read.table(paste(p, "_copynumber_segments.txt", sep=""), header=T)

# read two files written by segmenter; variables created are segLog and segBaf
load(paste(p, "_noSmear_segLog.RData", sep=""))
load(paste(p, "_noSmear_segBaf.RData", sep=""))

# get sample names
# MARY changed this line as it lost the first three samples
# sampNames <- intersect(colnames(logR_T[-c(1:3)]), colnames(segLog)[-c(1:5)])
sampNames <- intersect(colnames(logR_T[-c(1:2)]), colnames(segLog)[-c(1:2)])

getSegIdx <- function(idx, segs, pos) {
  seg=segs[idx,]
  subidx <- which(pos$Chr == seg$Chr)
  return(subidx[which(pos$Position[subidx] >= seg$Start & pos$Position[subidx] <= seg$End)])
}

segIdx <- lapply(1:nrow(cnSegs), getSegIdx, segs=cnSegs, pos=logR_T[,1:3])
seglen <- unlist(lapply(segIdx, length))
allSegIdx <- unlist(segIdx)

logR_T <- logR_T[allSegIdx, c("Chr","Position",sampNames)]
logR_N <- logR_N[allSegIdx, c("Chr","Position",sampNames)]
BAF_T <- BAF_T[allSegIdx, c("Chr","Position",sampNames)]
BAF_N <- BAF_N[allSegIdx, c("Chr","Position",sampNames)]

ascat.bc = ascat.loadData(logR_T, BAF_T, logR_N, BAF_N, chrs=1:22, read.from.disk=FALSE)
rm(list=c("BAF_T", "logR_T", "BAF_N", "logR_N"))
gc()

## Paste copynumber segmentation
## 3. Tumor_LogR_segmented: matrix of LogR segmented values
## 4. Tumor_BAF_segmented: list of BAF segmented values; each element in the list is a matrix containing the segmented values for one sample  (only for probes that are germline homozygous)

## logR
ascat.bc$segIndexes <- matrix(NA, nrow = nrow(cnSegs), ncol=2)
ascat.bc$Tumor_LogR_segmented <- matrix(NA, nc=length(sampNames), nr=nrow(ascat.bc$Tumor_LogR))
colnames(ascat.bc$Tumor_LogR_segmented) <- sampNames
rownames(ascat.bc$Tumor_LogR_segmented) <- rownames(ascat.bc$Tumor_LogR)
icount <- 1
for (i in 1:nrow(cnSegs)) {
  print(i)
  ascat.bc$Tumor_LogR_segmented[icount:(icount+length(segIdx[[i]])-1),] <-
    segLog[rep(i, length(segIdx[[i]])),colnames(ascat.bc$Tumor_LogR_segmented)]
  nextc <- icount + length(segIdx[[i]])
  ascat.bc$segIndexes[i,] = c(icount, nextc-1)
  icount <- nextc
}

## BAF
getGgBaf <- function(i, bafvals, seglen, ggidx, snames) {
  ##print(i)
  startidx <- sum(seglen[0:(i-1)])
  segidx <- intersect((startidx:(startidx+seglen[i]-1)), gg)
  if (length(segidx) > 0) {
    tmpv <- rep(bafvals[i], length(segidx))
    names(tmpv) <- snames[segidx]
    return(tmpv)
  }
  return(c())
}

ascat.bc$Tumor_BAF_segmented <- list()
for (i in 1:length(sampNames)) {
  s=sampNames[i]
  print(s)
  gg = which(ascat.bc$Germline_BAF[,s] >= 0.3 & ascat.bc$Germline_BAF[,s] <= 0.7)
  ascat.bc$Tumor_BAF_segmented[[i]] <- as.matrix(unlist(lapply(1:nrow(cnSegs), getGgBaf, bafvals=segBaf[,s], seglen=seglen, ggidx=gg,
                                                     snames=rownames(ascat.bc$SNPpos))))
}

save(ascat.bc, file=paste(p, "_fcn_ascat.bc.RData", sep=""))
##ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=0.55, allow100percent = T, doubled=doubled)
#ascat.output = ascat.runAscat(ascat.bc, gamma=0.55)
save(ascat.output, file=paste(p, "_fcn_ascat.output.RData", sep=""))

source("organize.ascat.segments.R")
snp.markers <- ascat.bc$SNPpos
ascat_segments <- organize.ascat.segments(ascat.output, snp.markers, baf=ascat.bc$Tumor_BAF,
                                                 logr=ascat.bc$Tumor_LogR, bafn=ascat.bc$Germline_BAF)
ascat_cont <- ascat.output$aberrantcellfraction
names(ascat_cont) <- colnames(ascat.output$nA)
ascat_ploidy <- ascat.output$ploidy
names(ascat_ploidy) <- colnames(ascat.output$nA)
ascat_psi <- ascat.output$psi
names(ascat_psi) <- colnames(ascat.output$nA)
ascat_goodness <- ascat.output$goodness
names(ascat_goodness) <- colnames(ascat.output$nA)

write.table(ascat_segments,file=paste(p, "_fcn_ascat_segments.txt", sep=""))

write.table(ascat_cont, file=paste(p, "_fcn_ascat_cont.txt", sep=""))

write.table(ascat_ploidy, file=paste(p, "_fcn_ascat_ploidy.txt", sep=""))

write.table(ascat.output$failedarrays, file=paste(p, "_failed_arrays.txt", sep=""), col.names=F, row.names=F, sep="\t", quote=F)

write.table(ascat_psi, file=paste(p, "_fcn_ascat_psi.txt", sep=""))

write.table(ascat_goodness, file=paste(p, "_fcn_ascat_goodness.txt", sep=""))
