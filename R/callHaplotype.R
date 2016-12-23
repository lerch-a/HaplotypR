
createContingencyTable <- function(inputFiles, dereplicated=F, inputFormat="fasta", outputDir=".", sampleNames, replicatNames, freqSplitPattern="_", haplotypeUIDFile=NULL){
  require(ShortRead)
  
  env <- environment()
  env$allHaplotypes <- DNAStringSet()
  if(!is.null(haplotypeUIDFile)){
    allHap <- readFasta(haplotypeUIDFile)
    env$allHaplotypes <- sread(allHap)
    names(env$allHaplotypes) <- id(allHap)
    rm(allHap)
  } 
  
  contingencyList <- lapply(seq_along(inputFiles), function(i){
    
    # load sequence reads
    if(inputFormat=="fasta")
      inputReads <- readFasta(inputFiles[i])
    else
      inputReads <- readFastq(inputFiles[i])
    
    # calculate amplicon frequencies
    if(!dereplicated){
      readFreq <- tables(sread(inputReads), n=NULL)$top
      haplotypes <- DNAStringSet(names(readFreq))
      names(readFreq) <- NULL
    }else{
      readFreq <- as.character(id(inputReads))
      readFreq <- as.integer(do.call(rbind, strsplit(readFreq, freqSplitPattern))[,2])
      haplotypes <- sread(inputReads)
    }
    
    # compare to all haplotypes, name haplotypes and add missing haplotypes to all haplotypes list
    idx <- match(haplotypes, env$allHaplotypes)
    inHapTab <- !is.na(idx)
    names(haplotypes)[inHapTab] <- names(env$allHaplotypes[idx[inHapTab]])
    names(haplotypes)[!inHapTab] <- paste("UID", seq_along(haplotypes[!inHapTab])+length(env$allHaplotypes), sep="")
    env$allHaplotypes <- append(env$allHaplotypes, haplotypes[!inHapTab])
    
    #
    #extractQualityPerHaplotype(haplotypes, inputReads, inputFiles[i])
    
    # return haplotyps per file
    lst <- list(UID=names(haplotypes), Frequency=readFreq, sampleName=sampleNames[i], replicatName=replicatNames[i])
    
    names(haplotypes) <- paste(names(haplotypes), readFreq, sep=freqSplitPattern)
    freqFilename <- file.path(outputDir, sprintf("%s_Rep%s_hapFreq.fa", sampleNames[i], replicatNames[i]))
    writeFasta(haplotypes, freqFilename)
    return(lst)
  })
  names(contingencyList) <- sprintf("%s_Rep%s", sampleNames, replicatNames)
  
  # return all haplotyps
  writeFasta(env$allHaplotypes, file.path(outputDir, "allHaplotypes.fa"))
  
  freq <- integer(length(env$allHaplotypes))
  names(freq) <- names(env$allHaplotypes)
  contingencyTable <- lapply(contingencyList, function(l){
    #idx <- match(l$UID, names(env$allHaplotypes))
    freq[as.character(l$UID)] <- l$Frequency
    return(freq)
  })
  contingencyTable <- do.call(cbind, contingencyTable)
  rownames(contingencyTable) <- names(env$allHaplotypes)
  colnames(contingencyTable) <- names(contingencyList) 
  
  return(contingencyTable)
}

callHaplotypePerSample <- function(x, sensitivity=1/1000, minCoverage=3L, minReplicate=2, reportBackground=T) {

  # # is sample replicated
  # hasReplicat <- dim(x)[2]>1
  
  # remove haplotypes without reads
  x <- x[rowSums(x>0) > 0,,drop=F]
  
  # selected and remove filtered haplotype
  idx <- rownames(x) %in% c("Chimera", "Singelton", "Indels", "Cut-Off_Sample", "Cut-Off_Size")
  background <- x[idx,, drop=F]
  x <- x[!idx,, drop=F]
  
  # check for noise haplotype
  minCov <- colSums(x)*sensitivity
  minCov[minCov<minCoverage] <- minCoverage
  noiseIdx <- rowSums(t(t(x)/minCov) >= 1) < minReplicate # only haplotypes present in both replicate
  lowCnt <- colSums(x[noiseIdx,,drop=F])
  x <- x[!noiseIdx,,drop=F]
  
  # add background to haplotyp counts
  if(reportBackground){
    if(dim(x)[1]>0)
      x <- rbind(Noise=lowCnt, background, x[order(rowSums(x)),,drop=F])
    else
      x <- t(data.frame(Noise=lowCnt, background))
  }

  
  return(x)
}

