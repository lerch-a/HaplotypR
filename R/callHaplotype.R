
createContingencyTable <- function(inputFiles, sampleNames, dereplicated=F, inputFormat="fasta", outputDir=".", replicatNames=NULL, freqSplitPattern="_", haplotypeUIDFile=NULL, progressReport=message){
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
    if(inputFormat=="fasta"){
      inputReads <- readFasta(inputFiles[i])
    }else{
      inputReads <- readFastq(inputFiles[i])
    }
    
    # remove reads with Ns
    idx <- grep("N", sread(inputReads), invert=T)
    inputReads <- inputReads[idx]
    rm(idx)
    
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
    freqFilename <- file.path(outputDir, sub(".fastq.gz$", "_hapFreq.fa", basename(inputFiles[i])))
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
  
  # if NA value replace with zero 
  x[is.na(x)] <- 0
  
  # remove low covarage sample
  cov <- colSums(x)
  idx <- cov>=minSampleCoverage
  if(all(!idx)){
    x[,!idx] <- NA
    x <- x[1,, drop=F]
    rownames(x) <- NA
    return(x)
  }else{
    x[,!idx] <- NA 
  }
  
  # sample replicated
  minReplicate <- min(dim(x)[2], minReplicate)
  
  # remove haplotypes without reads
  x <- x[rowSums(x>0, na.rm=T) > 0,,drop=F]
  
  # selected and remove filtered haplotype
  idx <- rownames(x) %in% defineBackground
  background <- x[idx,, drop=F]
  x <- x[!idx,, drop=F]
  
  # check for noise haplotype
  minCov <- colSums(x, na.rm=T)*detectability
  minCov[minCov<minHaplotypCoverage] <- minHaplotypCoverage
  if(dim(x)[1]>0){
    noiseIdx <- rowSums(t(t(x)/minCov) >= 1, na.rm=T) < minReplicate # only haplotypes present in minimum replicates
    lowCnt <- colSums(x[noiseIdx,,drop=F], na.rm=T)
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

createFinalHaplotypTable <- function(outputDir, sampleTable, markerTable, referenceSequence, snpList, postfix, 
                                     minHaplotypCoverage=3, minReplicate=2, 
                                     detectability=1/100, minSampleCoverage=300){
  # check args
  stopifnot(
    is.character(outputDir), length(outputDir) == 1, file.exists(outputDir),
    is.data.frame(sampleTable), all(c("MarkerID", "ReadFile", "SampleID") %in% colnames(sampleTable)),
    is.character(sampleTable$ReadFile), all(file.exists(sampleTable$ReadFile)),
    is.data.frame(markerTable), all(c("MarkerID") %in% colnames(markerTable)),
    is.list(snpList), all(sampleTable$MarkerID %in% names(snpList)),
    is.character(postfix), length(postfix) == 1,
    is.numeric(minHaplotypCoverage), length(minHaplotypCoverage) == 1,
    is.numeric(minReplicate), length(minReplicate) == 1,
    is.numeric(detectability), length(detectability) == 1,
    is.numeric(minSampleCoverage), length(minSampleCoverage) == 1
  )
  
  devMode <- getOption("HaplotypR.devel")
  if(is.null(devMode))
    devMode <- F
  
  outFreqDir <- file.path(outputDir, "frequencyFiles")
  dir.create(outFreqDir, recursive = T)
  res <- lapply(markerTable$MarkerID, function(marker){
    
    outFreqFiles <- file.path(outFreqDir, marker)
    dir.create(outFreqFiles)
    samTab <- sampleTable[sampleTable$MarkerID == marker,]
    snpSet <- as.data.frame(snpList[[marker]], stringsAsFactors = FALSE)
    snpSet$Pos <- as.integer(snpSet$Pos)
    snpSet$Alt <- NULL
    prefix <- sub(".fastq.gz$", "", basename(as.character(samTab$ReadFile)))
    
    
    # Create frequency files and count table
    tab <- createContingencyTable(inputFiles = as.character(samTab$ReadFile), sampleNames=as.character(samTab$SampleID), dereplicated=F,
                                  inputFormat="fastq", outputDir=outFreqFiles, replicatNames=NULL, haplotypeUIDFile=NULL)
    if(devMode)
      write.table(tab, file=file.path(outputDir, sprintf("contingencyTable_%s%s.txt", marker, postfix)), sep="\t")
    fnAllSeq <- file.path(outFreqFiles, sprintf("allSequences_%s%s.fasta", marker, postfix))
    file.rename(file.path(outFreqFiles, "allHaplotypes.fa"), fnAllSeq)
    frqfile <- file.path(outFreqFiles, paste(prefix, "_hapFreq.fa", sep=""))
    
    # run cluster with Rswarm package
    outCluster <- file.path(outputDir, "cluster", marker)
    dir.create(outCluster, recursive=T)
    clusterFilenames <- clusterReads(frqfile, outCluster, prefix)
    
    # check for chimeras with Rvsearch package
    chimeraFilenames <- checkChimeras(clusterFilenames[,"RepresentativeFile"], method="vsearch")
    
    # create an overview table
    overviewHap <- createHaplotypOverviewTable(allHaplotypesFilenames=fnAllSeq, clusterFilenames=clusterFilenames, 
                                               chimeraFilenames=chimeraFilenames, referenceSequence=referenceSequence[marker],
                                               snpSet=snpSet)
    
    ## create final haplotype
    overviewHap$FinalHaplotype <- NA_character_
    overviewHap[overviewHap$representatives, "FinalHaplotype"] <- marker
    # Set singelton
    overviewHap[overviewHap$singelton, "FinalHaplotype"] <- "Singelton"
    # Set Indel
    overviewHap[!is.na(overviewHap$indels) & overviewHap$indels & overviewHap$representatives, "FinalHaplotype"] <- "Indels"
    # Set chimera
    overviewHap[!is.na(overviewHap$chimeraScore) & is.na(overviewHap$nonchimeraScore), "FinalHaplotype"] <- "Chimera"
    # Cluster identical SNPs pattern
    idx <- overviewHap$FinalHaplotype %in% marker
    snps <- unique(na.omit(overviewHap$snps[idx]))
    if(length(snps) > 0)
      names(snps) <- paste(marker, seq_along(snps), sep="-")
    overviewHap$FinalHaplotype[idx] <- names(snps)[match(overviewHap$snps[idx], snps)]
    overviewHap <- as.data.frame(overviewHap)
    if(devMode)
      write.table(overviewHap, file=file.path(outputDir, sprintf("HaplotypeOverviewTable_%s%s.txt", marker, postfix)), sep="\t")
    
    # getHaplotype sequence
    if(!is.null(referenceSequence)){
      hapSeq <- lapply(snps, function(sp){
        replaceLetterAt(referenceSequence[[marker]], at=snpSet$Pos, letter=sp)
      })
      hapSeq <- DNAStringSet(hapSeq)
      writeFasta(hapSeq, file.path(outputDir, file=sprintf("%s_HaplotypeSeq%s.fasta", marker, postfix)))
    }
    repfile <- clusterFilenames[,"RepresentativeFile"]
    
    hab <- rep(0, dim(overviewHap)[1])
    names(hab) <- rownames(overviewHap)
    
    haplotypesSample <- lapply(seq_along(repfile), function(i){
      sr1 <- readFasta(repfile[i])
      vec <- do.call(rbind, strsplit(as.character(id(sr1)), "_"))
      clusterResp <- vec[,1]
      clusterSize <- as.integer(vec[,2])
      hab[clusterResp] <- clusterSize
      hab
    })
    haplotypesSample <- do.call(cbind, haplotypesSample)
    haplotypesSample <- haplotypesSample[rowSums(haplotypesSample)>0, , drop = FALSE]
    colnames(haplotypesSample) <- sub(".representatives.fasta", "", basename(repfile))
    
    overviewHap <- overviewHap[rownames(haplotypesSample),]
    
    # set final haplotype names
    rownames(haplotypesSample) <- overviewHap[rownames(haplotypesSample), "FinalHaplotype"]
    if(devMode)
      write.table(cbind(HaplotypNames=rownames(haplotypesSample),haplotypesSample), 
                  file=file.path(outputDir, sprintf("rawHaplotypeTable_%s%s.txt", marker, postfix)), sep="\t", row.names=F, col.names=T)
    
    haplotypesSample <- reclusterHaplotypesTable(haplotypesSample)
    if(devMode)
      write.table(cbind(HaplotypNames=rownames(haplotypesSample),haplotypesSample), 
                  file=file.path(outputDir, sprintf("reclusteredHaplotypeTable_%s%s.txt", marker, postfix)), sep="\t", row.names=F, col.names=T)
    
    
    # Apply cut-off haplotype only in 1 sample
    haplotypesSample <- callHaplotype(haplotypesSample, minHaplotypCoverage=minHaplotypCoverage, minReplicate=minReplicate, 
                                      detectability=detectability, minSampleCoverage=1)
    
    if(devMode) 
      write.table(cbind(HaplotypNames=rownames(haplotypesSample), haplotypesSample), 
                  file=file.path(outputDir, sprintf("finalHaplotypeTable_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt",
                                                    minCov, minCovSample, minReplicate, detectability, marker, postfix)),
                  sep="\t", row.names=F, col.names=T)
    
    
    # check replicates
    idx <- split(1:dim(haplotypesSample)[2], samTab$SampleID)
    markerRes <- lapply(idx, function(i){
      reads <- callHaplotype(haplotypesSample[,i, drop=F], minHaplotypCoverage=minHaplotypCoverage, 
                             minReplicate=minReplicate, detectability=detectability, minSampleCoverage=minSampleCoverage, 
                             reportBackground=T)
      tab <- cbind(samTab[rep(i,each=dim(tab)[1]), c("SampleID","SampleName","MarkerID")], 
                   Haplotype=rownames(tab), Reads=as.integer(tab), FlagChimera=F)
      colnames(tab) <- c("SampleID","SampleName","MarkerID","Haplotype","Reads")
      rownames(tab) <- NULL
      
      #check individual samples for chimera
      do.call(rbind, lapply(split(tab, tab$SampleID), function(tt){
        chim <- NULL
        hIdx <- grep(marker, tt$Haplotype)
        if(length(hIdx)>2){
          chim <- flagChimera(tt[hIdx,], overviewHap)
        }
        tt$FlagChimera <- tt$Haplotype %in% chim
        return(tt)
      }))
      return(tab)
    })
    markerRes <- do.call(rbind.data.frame, markerRes)
    rownames(markerRes) <- NULL
    markerResFn <- file.path(outputDir, sprintf("finalHaplotypeList_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt", 
                                                minHaplotypCoverage, minSampleCoverage, minReplicate, detectability, marker, postfix))
    write.table(markerRes, file=markerResFn, sep="\t", row.names=F, col.names=T)
    return(markerRes)
  })
  names(res) <- markerTable$MarkerID
  return(res)
}
