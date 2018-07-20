organize.ascat.segments <- function(ascat.output, markers, baf=NULL, logr=NULL, bafn=NULL, dev05=0.2){
  ## markers = matrix with two columns with SNP positions
  samp.names <- colnames(ascat.output$nA)
  failed.names <- ascat.output$failedarrays
  if(length(failed.names) >= 1){
    failed.locations <- which(is.na(ascat.output$segments))
    new.samp.names <- 1:length(ascat.output$segments)
    new.samp.names[-c(failed.locations)] <- samp.names
    new.samp.names[failed.locations] <- failed.names
    samp.names <- new.samp.names
  }
  out.seg <- matrix(0,0,ncol=9)
  colnames(out.seg) <- c("SampleID", "Chr", "Start", "End", "nProbes", "mBAF", "logR", "nA", "nB")
  if(length(ascat.output$segments) != length(samp.names)){stop}
  for(i in 1:length(ascat.output$segments)){
    tmp.out <- out.seg[0,]
    sample.segs <- ascat.output$segments[[i]]
    print.noquote(paste("Processing sample",samp.names[i]))
    if(class(sample.segs) != "matrix"){ next }
    if(class(sample.segs) == "matrix"){
      isErr <- 0
      for(j in 1:nrow(sample.segs)){
        if (isErr) {next}
        tmp <- markers[sample.segs[j,1]:sample.segs[j,2],,drop=F]
        tmpbaf <- tmplogr <- rep(NA,nrow(tmp))
        if (!is.null(baf)) {
          tmpbaf <- baf[sample.segs[j,1]:sample.segs[j,2],samp.names[i]]
          if (!is.null(bafn)) {
            tmpbafn <- bafn[sample.segs[j,1]:sample.segs[j,2],samp.names[i]]
            tmpbaf[which(tmpbafn <= (0.5 - dev05) | tmpbafn >= (0.5 + dev05))] <- NA
          } else {
            tmpbaf[which(tmpbaf <= 0.15 | tmpbaf >= 0.85)] <- NA
          }
          tmpbaf <- 0.5 + abs(0.5 - tmpbaf)
        }
        if (!is.null(logr)) {
          tmplogr <- logr[sample.segs[j,1]:sample.segs[j,2],samp.names[i]]
        }
        if(!all(tmp[,1]== tmp[1,1])){
          if(length(unique(tmp[,1])) >2) {
            cat("Error in segmentation, more than two chromosomes in one segment. Skipping\n")
            isErr <- 1
            next
          }

          cat("Error in segmentation, two chromosomes in one segment. Splitting into parts.\n") # ASCAT joins chromosome regions smaller than 200 probes to the next one
          tmp1 <- tmp[tmp[,1] == tmp[1,1],]
          tmpbaf1 <- tmpbaf[tmp[,1] == tmp[1,1]]
          tmplogr1 <- tmplogr[tmp[,1] == tmp[1,1]]
          tmp2 <- tmp[tmp[,1] == tmp[nrow(tmp),1],]
          tmpbaf2 <- tmpbaf[tmp[,1] == tmp[nrow(tmp),1]]
          tmplogr2 <- tmplogr[tmp[,1] == tmp[nrow(tmp),1]]
          chr1 <- as.character(tmp1[1,1])
          seg.start1 <- as.character(tmp1[1,2])
          seg.end1 <- as.character(tmp1[nrow(tmp1),2])
          nProbes1 <- nrow(tmp1)
          tmp.out <- rbind(tmp.out, c(samp.names[i], chr1, seg.start1, seg.end1, nProbes1,
                                      mean(tmpbaf1,na.rm=T), mean(tmplogr1, na.rm=T), sample.segs[j,3:4]))

          chr2 <- as.character(tmp2[1,1])
          seg.start2 <- as.character(tmp2[1,2])
          seg.end2 <- as.character(tmp2[nrow(tmp2),2])
          nProbes2 <- nrow(tmp2)
          tmp.out <- rbind(tmp.out, c(samp.names[i], chr2, seg.start2, seg.end2, nProbes2,
                                      mean(tmpbaf2, na.rm=T), mean(tmplogr2, na.rm=T), sample.segs[j,3:4]))
          
        }
        if(all(tmp[,1]== tmp[1,1])){
          chr <- as.character(tmp[1,1])
          seg.start <- as.character(tmp[1,2])
          seg.end <- as.character(tmp[nrow(tmp),2])
          nProbes <- nrow(tmp)
          tmp.out <- rbind(tmp.out, c(samp.names[i], chr, seg.start, seg.end, nProbes,
                                      mean(tmpbaf, na.rm=T), mean(tmplogr, na.rm=T), sample.segs[j,3:4]))
        }
      }
      if (!isErr) { out.seg <- rbind(out.seg, tmp.out) }
      gc()
    }
  }
  out.seg[out.seg[,2] == "X",2] <- 23
  out.seg <- as.data.frame(out.seg)
  out.seg[,2:9] <- apply(out.seg[,2:9], 2, function(x){as.numeric(as.character(x))})
  out.seg[,1] <- as.character(out.seg[,1])
  return(out.seg)
}
