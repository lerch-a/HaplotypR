
createContingencyTable <- function(inputFiles, dereplicated=F, inputFormat="fasta", outputDir=".", sampleNames, replicatNames=NULL, freqSplitPattern="_", haplotypeUIDFile=NULL, progressReport=message){
  require(ShortRead)
  require(Biostrings)
  
  if(is.null(replicatNames))
    replicatNames <- rep("", length(sampleNames))

  env <- environment()
  env$allHaplotypes <- DNAStringSet()
  if(!is.null(haplotypeUIDFile)){
    allHap <- readFasta(haplotypeUIDFile)
    env$allHaplotypes <- sread(allHap)
    names(env$allHaplotypes) <- id(allHap)
    rm(allHap)
  } 

  contingencyList <- lapply(seq_along(inputFiles), function(i){
    
    # check and set progress report function
    if(!is.function(progressReport))
      progressReport <- message
    msg <- paste("Processing file", basename(inputFiles[i]), "...", sep=" ")
    progressReport(detail=msg, value=i)
    
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
    idx <- Biostrings::match(haplotypes, env$allHaplotypes)
    inHapTab <- !is.na(idx)
    names(haplotypes)[inHapTab] <- names(env$allHaplotypes[idx[inHapTab]])
    names(haplotypes)[!inHapTab] <- paste("UID", seq_along(haplotypes[!inHapTab])+length(env$allHaplotypes), sep="")
    env$allHaplotypes <- append(env$allHaplotypes, haplotypes[!inHapTab])
    
    # return haplotyps per file
    lst <- list(UID=names(haplotypes), Frequency=readFreq, sampleName=sampleNames[i], replicatName=replicatNames[i])
    
    names(haplotypes) <- paste(names(haplotypes), readFreq, sep=freqSplitPattern)
    #freqFilename <- file.path(outputDir, sprintf("%s%s_hapFreq.fa", sampleNames[i], replicatNames[i]))
    freqFilename <- file.path(outputDir, sub(".fastq.gz","_hapFreq.fa", basename(inputFiles[i])))
    writeFasta(haplotypes, freqFilename)
    return(lst)
  })
  names(contingencyList) <- sprintf("%s%s", sampleNames, replicatNames)
  
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

callHaplotype <- function(x, detectability=1/100, minHaplotypCoverage=3, minReplicate=2, minSampleCoverage=300, reportBackground=T, defineBackground=NULL, ...) {

  # check minHaplotypCoverage argument
  if(minHaplotypCoverage < 3){
    arguments <- list(...)
    if(any(is.null(arguments$overwriteMinCoverage), arguments$overwriteMinCoverage==F)){
      stop("The minimum read coverage per Haplotype (= minHaplotypCoverage) must be at least 3. 
           To overwrite this default minHaplotypCoverage setting see details section of manual.")
    }
  }
  
  # check defineBackground argument
  if(is.null(defineBackground))
    defineBackground <- c("Chimera", "Singelton", "Indels", "Cut-Off_Sample", "Cut-Off_Size")
  
 	# remove low covarage sample
  cov <- colSums(x)
  idx <- cov>minSampleCoverage
  if(all(!idx)){
  	x[,idx] <- NA
  	rownames(x) <- NA
  	return(x[,idx])
  }else{
  	x[,idx] <- NA 
  }

  # sample replicated
  minReplicate <- min(dim(x)[2], minReplicate)
  
  # remove haplotypes without reads
  x <- x[rowSums(x>0) > 0,,drop=F]
  
  # selected and remove filtered haplotype
  idx <- rownames(x) %in% defineBackground
  background <- x[idx,, drop=F]
  x <- x[!idx,, drop=F]
  
  # check for noise haplotype
  minCov <- colSums(x)*detectability
  minCov[minCov<minHaplotypCoverage] <- minHaplotypCoverage
  if(dim(x)[1]>0){
  	noiseIdx <- rowSums(t(t(x)/minCov) >= 1) < minReplicate # only haplotypes present in minimum replicates
  	lowCnt <- colSums(x[noiseIdx,,drop=F])
  	x <- x[!noiseIdx,,drop=F]
  }else{
  	lowCnt <- 0
  }
  
  # add background to haplotyp counts
  if(reportBackground){
    if(dim(x)[1]>0)
      x <- rbind(Noise=lowCnt, background, x[order(rowSums(x)),,drop=F])
    else
    	x <- rbind(Noise=lowCnt, background)
  }

  return(x)
}

