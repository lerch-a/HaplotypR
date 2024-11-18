
rmNreads <- function(x) {
  numN <- alphabetFrequency(sread(x), baseOnly=T)[,"other"]
  return(x[numN==0])
}


createFinalHaplotypTableDADA2 <- function(outputDir, sampleTable, markerTable, referenceSequence=NULL, snpList=NULL, postfix="", 
                                     minHaplotypCoverage=3, minReplicate=2, 
                                     detectability=1/100, minSampleCoverage=300,
                                     minSeqLength=5, maxSeqLength=9999, filterIndel=T,
                                     tempdir=tempdir(),
                                     multithread=FALSE, pool="pseudo", OMEGA_A=1e-120){
  require(dada2)
  
  # check if input file exists
  sampleTable <- sampleTable[!is.na(sampleTable$ReadFile),]
  sampleTable <- sampleTable[file.exists(sampleTable$ReadFile),]
  
  # check args
  stopifnot(
    is.character(outputDir), length(outputDir) == 1, file.exists(outputDir),
    is.data.frame(sampleTable), all(c("MarkerID", "ReadFile", "SampleID", "SampleName","numRead") %in% colnames(sampleTable)),
    is.character(sampleTable$ReadFile),
    is.data.frame(markerTable), all(c("MarkerID") %in% colnames(markerTable)),
    is.character(postfix), length(postfix) == 1,
    is.numeric(minHaplotypCoverage), length(minHaplotypCoverage) == 1,
    is.numeric(minReplicate), length(minReplicate) == 1,
    is.numeric(detectability), length(detectability) == 1,
    is.numeric(minSampleCoverage), length(minSampleCoverage) == 1
  )
  if(is.null(referenceSequence)){
    referenceSequence <- markerTable$ReferenceSequence
    names(referenceSequence) <- markerTable$MarkerID
  }
  
  devMode <- getOption("HaplotypR.devel")
  if(is.null(devMode))
    devMode <- F
  
  if(is.null(referenceSequence)){
    referenceSequence <- DNAStringSet(markerTable$ReferenceSequence)
    names(referenceSequence) <- markerTable$MarkerID
  }

  # filter read file
  filtDir <- file.path(tempdir,"filtered")
  dir.create(filtDir)
  newFile <- file.path(filtDir, basename(sampleTable$ReadFile))
  numReads <- unlist(lapply(seq_along(sampleTable$ReadFile), function(ii){
    sr <- readFastq(sampleTable$ReadFile[ii])
    sr <- sr[minSeqLength<width(sr) & width(sr)<maxSeqLength]
    sr <- rmNreads(sr)
    #length(sr)
    writeFastq(sr, newFile[[ii]], compress = T)
  }))
  sampleTable$ReadFile <- newFile
  sampleTable$numReads <- numReads
  sampleTable <- sampleTable[sampleTable$numReads>=3,]
  #write.table(sampleTable, "demultiplexMarkerSummary_MinION_filt.txt", sep="\t", row.names=F)
  
  # out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
  #                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
  #                      compress=TRUE, multithread=TRUE)
  
  # run dada2
  err <- learnErrors(sampleTable$ReadFile, multithread=multithread)
  # plotErrors(err, nominalQ=TRUE)
  dadas <- dada(sampleTable$ReadFile, err=err, multithread=multithread, pool=pool, OMEGA_A=OMEGA_A)
  seqtab <- makeSequenceTable(dadas)
  # remove chimera
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=TRUE)
  #sum(seqtab.nochim)/sum(seqtab)
  seqtabLst <- split(as.data.frame(seqtab.nochim), sampleTable$MarkerID)
  file.remove(newFile)
  
  selMarker <- unique(sampleTable$MarkerID)
  resultsLst <- lapply(selMarker, function(nm){
    # nm <- "ama1_D3"
    # message(nm)
    # get read counts and rename sample
    markTab <- seqtabLst[[nm]]
    markTab <- markTab[,colSums(markTab)>0, drop=F]
    rownames(markTab) <- sampleTable$SampleID[match(rownames(markTab), basename(sampleTable$ReadFile))]
    # get haplotype sequence
    haplotyp <- DNAStringSet(colnames(markTab))
    
    # check sample cutoff
    totReads <- rowSums(markTab, na.rm=T)
    # summary(totReads)
    Noise1 <- markTab
    minCovSam <- totReads>minSampleCoverage
    markTab[!minCovSam,] <- 0 # apply threshold!
    Noise1[minCovSam,] <- 0 # apply threshold!
    # check minority clone cutoff
    Noise2 <- markTab
    detectionLimit <- (markTab/totReads)>=detectability
    minCovHap <- markTab>minHaplotypCoverage
    markTab[!(detectionLimit & minCovHap)] <- 0 # apply threshold!
    Noise2[detectionLimit | minCovHap] <- 0 
    keepHaplotype <- colSums(markTab)>0
    markTab <- markTab[,keepHaplotype, drop=F]
    haplotyp <- haplotyp[keepHaplotype]
    # check for length polymorphism
    if(filterIndel){
      mLength <- (width(referenceSequence[nm])-1)<width(haplotyp) & 
        width(haplotyp)<(width(referenceSequence[nm])+1) 
      Noise3 <- markTab[,!mLength, drop=F]
      markTab <- markTab[,mLength, drop=F]
      haplotyp <- haplotyp[mLength]
    }
    # rename and save haplotyp sequence
    if(length(haplotyp)>0){
      names(haplotyp) <- paste(nm, "_Id", 1:length(haplotyp), sep="")
      colnames(markTab) <- names(haplotyp)
    }
    # save number of censored reads
    markTab <- cbind(markTab, Censored=totReads-rowSums(markTab, na.rm=T))
    writeFasta(haplotyp, file=file.path(outputDir, sprintf("HaplotypeList_%s.fasta", nm)))
    return(list(markTab, haplotyp))
  })
  names(resultsLst) <- selMarker

  file.remove(newFile[file.exists(newFile)])
  file.remove(filtDir)
  finalHaplotyopList <- writeHaplotypList(resultsLst)
  #write.csv(finalHaplotyopList, file=file.path(outputDir, "finalHaplotypList_MinION_MH.csv"), row.names=F)
  
  ## Start TODO
  #browser()

  # set final haplotype names
  # rownames(haplotypesSample) <- overviewHap[rownames(haplotypesSample), "FinalHaplotype"]
  # if(devMode)
  #   write.table(cbind(HaplotypNames=rownames(haplotypesSample),haplotypesSample), 
  #               file=file.path(outputDir, sprintf("rawHaplotypeTable_%s%s.txt", marker, postfix)), sep="\t", row.names=F, col.names=T)
  
  # haplotypesSample <- reclusterHaplotypesTable(haplotypesSample)
  # if(devMode)
  #   write.table(cbind(HaplotypNames=rownames(haplotypesSample),haplotypesSample), 
  #               file=file.path(outputDir, sprintf("reclusteredHaplotypeTable_%s%s.txt", marker, postfix)), sep="\t", row.names=F, col.names=T)
  
  # # Apply cut-off haplotype only in 1 sample
  # haplotypesSample <- callHaplotype(haplotypesSample, minHaplotypCoverage=minHaplotypCoverage, minReplicate=minReplicate, 
  #                                   detectability=detectability, minSampleCoverage=1)
  
  # if(devMode) 
  #   write.table(cbind(HaplotypNames=rownames(haplotypesSample), haplotypesSample), 
  #               file=file.path(outputDir, sprintf("finalHaplotypeTable_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt",
  #                                                 minHaplotypCoverage, minSampleCoverage, minReplicate, detectability, marker, postfix)),
  #               sep="\t", row.names=F, col.names=T)
  
  
  # # check replicates
  # idx <- split(1:dim(haplotypesSample)[2], samTab$SampleID)
  # markerRes <- lapply(idx, function(i){
  #   tab <- callHaplotype(haplotypesSample[,i, drop=F], minHaplotypCoverage=minHaplotypCoverage, 
  #                        minReplicate=minReplicate, detectability=detectability, minSampleCoverage=minSampleCoverage, 
  #                        reportBackground=T)
  #   tab <- cbind(samTab[rep(i,each=dim(tab)[1]), c("SampleID","SampleName","MarkerID")], 
  #                Haplotype=rownames(tab), Reads=as.integer(tab), FlagChimera=F)
  #   colnames(tab) <- c("SampleID","SampleName","MarkerID","Haplotype","Reads")
  #   rownames(tab) <- NULL
  #   
  #   #check individual samples for chimera
  #   do.call(rbind, lapply(split(tab, tab$SampleID), function(tt){
  #     chim <- NULL
  #     hIdx <- grep(marker, tt$Haplotype)
  #     if(length(hIdx)>2){
  #       chim <- flagChimera(tt[hIdx,], overviewHap)
  #     }
  #     tt$FlagChimera <- tt$Haplotype %in% chim
  #     return(tt)
  #   }))
  #   return(tab)
  # })
  # markerRes <- do.call(rbind.data.frame, markerRes)
  # rownames(markerRes) <- NULL
  # markerResFn <- file.path(outputDir, sprintf("finalHaplotypeList_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt", 
  #                                             minHaplotypCoverage, minSampleCoverage, minReplicate, detectability, marker, postfix))
  # write.table(markerRes, file=markerResFn, sep="\t", row.names=F, col.names=T)
  # return(markerRes)
  ## End TODO
  return(finalHaplotyopList)
}

