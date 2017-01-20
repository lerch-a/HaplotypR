
calculateMismatchFrequencies <- function(fastqFiles, referenceSequence, excludeSNPList=NULL, minCoverage=50L, progressReport=message){
  
  seqErrC <- lapply(seq_along(fastqFiles), function(i){
    # check and set progress report function
    if(!is.function(progressReport))
      progressReport <- message
    msg <- paste("Processing file", basename(fastqFiles[i]), "...", sep=" ")
    progressReport(detail=msg, value=i)
    
    sr1 <- readFastq(fastqFiles[i])
    
    if(length(sr1)>=minCoverage){
      aln1 <- pairwiseAlignment(sread(sr1), referenceSequence, type="global")
      # idx <- (lengths(deletion(aln1)) + lengths(insertion(aln1))) == 0
      # table(idx)
      
      seqErr <- tapply(mismatchSummary(aln1)$subject$Count, mismatchSummary(aln1)$subject$SubjectPosition, sum)
      
      # exclude SNP from SNPList
      if(!is.null(excludeSNPList)){
        idx <- paste(mismatchSummary(aln1)$subject$Pattern, mismatchSummary(aln1)$subject$SubjectPosition) %in% 
          paste(excludeSNPList[,"Alt"], excludeSNPList[,"Pos"])
      }else{
        idx <- rep(F, length(mismatchSummary(aln1)$subject$Count))
      }
      
      seqErr <- tapply(mismatchSummary(aln1)$subject$Count[!idx], mismatchSummary(aln1)$subject$SubjectPosition[!idx], sum)
      seqErrF <- rep(0, width(referenceSequence))
      seqErrF[as.integer(names(seqErr))] <- seqErr
      
      return(cbind(MisMatch=seqErrF, Coverage=as.vector(coverage(aln1))))
    } else return(NULL)
  })
  names(seqErrC) <- fastqFiles
  seqErrC

}

callGenotype <- function(mismatchRateTable, minMismatchRate=0.5, minReplicate=2){
  potSNP <- rowSums(mismatchRateTable>minMismatchRate)>=minReplicate
  potSNP <- seq_along(potSNP)[potSNP]
  return(potSNP)
}
