
calculateMismatchFrequencies <- function(fastqFiles, referenceSequence, method=c("compareDNAString","pairwiseAlignment"), excludeSNPList=NULL, minCoverage=50L, progressReport=message){
  
  seqErrC <- lapply(seq_along(fastqFiles), function(i){
    # check and set progress report function
    if(!is.function(progressReport))
      progressReport <- message
    msg <- paste("Processing file", basename(fastqFiles[i]), "...", sep=" ")
    progressReport(detail=msg, value=i)
    
    sr1 <- readFastq(fastqFiles[i])
    
    if(length(sr1)>=minCoverage){
    	if(method == "pairwiseAlignment"){
    		aln1 <- pairwiseAlignment(sread(sr1), referenceSequence, type="global")
    		mismatchSummary <- mismatchSummary(aln1)$subject    		
    	} else
    		mismatchSummary <- compareDNAString(sread(sr1), referenceSequence)
      # idx <- (lengths(deletion(aln1)) + lengths(insertion(aln1))) == 0
      # table(idx)
      
      seqErr <- tapply(mismatchSummary$Count, mismatchSummary$SubjectPosition, sum)
      
      # exclude SNP from SNPList
      if(!is.null(excludeSNPList)){
        idx <- paste(mismatchSummary$Pattern, mismatchSummary$SubjectPosition) %in% 
          paste(excludeSNPList[,"Alt"], excludeSNPList[,"Pos"])
      }else{
        idx <- rep(F, length(mismatchSummary$Count))
      }
      
      seqErr <- tapply(mismatchSummary$Count[!idx], mismatchSummary$SubjectPosition[!idx], sum)
      seqErrF <- rep(0, width(referenceSequence))
      seqErrF[as.integer(names(seqErr))] <- seqErr
      
      return(cbind(MisMatch=seqErrF, Coverage=as.vector(length(sr1))))
    } else return(NULL)
  })
  names(seqErrC) <- fastqFiles
  seqErrC

}

callGenotype <- function(mismatchRateTable, minMismatchRate=0.5, minReplicate=2){
  potSNP <- rowSums(mismatchRateTable>minMismatchRate, na.rm=T)>=minReplicate
  potSNP <- seq_along(potSNP)[potSNP]
  return(potSNP)
}
