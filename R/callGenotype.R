
calculateMismatchFrequencies <- function(fastqFiles, referenceSequence, method=c("pairwiseAlignment", "compareDNAString"), excludeSNPList=NULL, minCoverage=50L, progressReport=message){

  method <- match.arg(method)
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
    		mismatchSummary <- pwalign::mismatchSummary(aln1)$subject    		
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


createSNPsList <-  function(outputDir, sampleTable, markerTable, refSeq, postfix="",
                                          method=c("pairwiseAlignment", "compareDNAString"), 
                                          minMMrate=0.5, minOccGen=2, minCoverage=50L,
                                          excludeSNPList=NULL, progressReport=message){
  # fix for historical inconsistancy
  if(!"ReadFile" %in% colnames(sampleTable)){
    sampleTable$ReadFile <- sampleTable$FileR1
  }
  
  # check if input file exists
  sampleTable <- sampleTable[!is.na(sampleTable$ReadFile),]
  sampleTable <- sampleTable[file.exists(sampleTable$ReadFile),]
  
  # Calculate mismatch rate
  seqErrLst <- lapply(markerTab$MarkerID, function(marker){
    seqErrorLst <- calculateMismatchFrequencies(as.character(sampleTable[sampleTable$MarkerID == marker, "ReadFile"]), 
                                              refSeq[marker], 
                                              method ="pairwiseAlignment", # c("pairwiseAlignment","compareDNAString"), 
                                              minCoverage=100L)
    names(seqErrorLst) <- sampleTable[sampleTable$MarkerID == marker, "SampleID"]
    seqErr <- do.call(cbind, lapply(seqErrorLst, function(l){
      l[,"MisMatch"]/l[,"Coverage"]
    }))
    write.table(seqErr, file.path(outputDir, sprintf("mismatchRate_rate_%s%s.txt", marker, postfix)), sep="\t", row.names=F)
    return(seqErr)
  })
  names(seqErrLst) <- markerTab$MarkerID
  
  # Call SNPs
  snpLst <- lapply(markerTab$MarkerID, function(marker){
    seqErr <- seqErrLst[[marker]]
    potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
    snpRef <- unlist(lapply(potSNP, function(snp){
      as.character(subseq(refSeq[marker], start=snp, width=1))
    }))
    snps <- data.frame(Chr=marker, Pos=potSNP, Ref=snpRef, Alt="N", stringsAsFactors=F)
    write.table(snps, file=file.path(outputDir, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt",
                                                    minMMrate*100, minOccGen, marker, postfix)),
                row.names=F, col.names=T, sep="\t", quote=F)
    return(snps)
  })
  names(snpLst) <- markerTab$MarkerID

  # Plot mismatch rate and SNP calls
  invisible(lapply(markerTab$MarkerID, function(marker){
    png(file.path(outputDir, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png",
                                     minMMrate*100, minOccGen, marker, postfix)),
        width=1500 , height=600)
    matplot(seqErrLst[[marker]], type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1),
            ylab="Mismatch Rate", xlab="Base Position", main=marker, cex.axis=2, cex.lab=2)
    abline(v=snpLst[[marker]][,"Pos"], lty=2, col="grey")
    abline(h=minMMrate, lty=1, col="red")
    dev.off()
  }))
  return(snpLst)
}
