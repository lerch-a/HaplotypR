
calculateMismatchFrequencies <- function(fastqFiles, referenceSequence, excludeSNPList=NULL, minCoverage=50L, ...){
  
  seqErrC <- lapply(seq_along(fastqFiles), function(i){
    message("Processing file", fastqFiles[i], "...", sep=" ")
    sr1 <- readFastq(fastqFiles[i])
    
    if(length(sr1)>=minCoverage){
      aln1 <- pairwiseAlignment(sread(sr1), referenceSequence, type="global")
      # idx <- (lengths(deletion(aln1)) + lengths(insertion(aln1))) == 0
      # table(idx)
      
      seqErr <- tapply(mismatchSummary(aln1)$subject$Count, mismatchSummary(aln1)$subject$SubjectPosition, sum)
      
      # exclude known FC27 snp
      idx <- logical()
      if(!is.null(excludeSNPList)){
        idx <- paste(mismatchSummary(aln1)$subject$Pattern, mismatchSummary(aln1)$subject$SubjectPosition) %in% 
          paste(mismatchSummary(aln)$subject$Pattern, mismatchSummary(aln)$subject$SubjectPosition)
      }
      
      seqErr <- tapply(mismatchSummary(aln1)$subject$Count[!idx], mismatchSummary(aln1)$subject$SubjectPosition[!idx], sum)
      seqErrF <- rep(0, width(referenceSequence))
      seqErrF[as.integer(names(seqErr))] <- seqErr
      
      return(cbind(MisMatch=seqErrF, Coverage=as.vector(coverage(aln1))))
    } else return(NULL)
  })
  names(seqErrC) <- fastqFiles

}

