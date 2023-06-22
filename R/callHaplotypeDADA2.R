
writeHaplotypList <- function(resLst){
  do.call(rbind, lapply(names(resLst), function(nm){
    message(nm)
    do.call(rbind, lapply(rownames(resLst[[nm]][[1]]), function(id){
      # nm <- names(resLst)[1]
      # id <- rownames(resLst[[nm]][[1]])[1]
      hap <- unlist(resLst[[nm]][[1]][id,,drop=F])
      samId <- strsplit(id, "_BC")[[1]][1]
      df <- data.frame(SampleID=samId, MarkerID=nm, 
                       Haplotype=names(hap),
                       Reads=hap)
      if(length(resLst[[nm]])>=3)
        df$Strain=resLst[[nm]][[3]][names(hap)]
      else
        df$Strain=NA
      if(length(resLst[[nm]])==4)
        df$Mutations=resLst[[nm]][[4]][names(hap)]
      else
        df$Mutations=NA
      # remove haplotyops with 0 reads
      df <- df[df$Reads>0,,drop=F]
      return(df)
    }))
  }))
}

createFinalHaplotypTableDADA2 <- function(outputDir, sampleTable, markerTable, referenceSequence=NULL, snpList=NULL, postfix, 
                                     minHaplotypCoverage=3, minReplicate=2, 
                                     detectability=1/100, minSampleCoverage=300,
                                     filterIndel=T){
  required(dada2)
  
  # check args
  stopifnot(
    is.character(outputDir), length(outputDir) == 1, file.exists(outputDir),
    is.data.frame(sampleTable), all(c("MarkerID", "ReadFile", "SampleID", "SampleName") %in% colnames(sampleTable)),
    is.character(sampleTable$ReadFile), all(file.exists(sampleTable$ReadFile)),
    is.data.frame(markerTable), all(c("MarkerID") %in% colnames(markerTable)),
    is.list(snpList),
    is.character(postfix), length(postfix) == 1,
    is.numeric(minHaplotypCoverage), length(minHaplotypCoverage) == 1,
    is.numeric(minReplicate), length(minReplicate) == 1,
    is.numeric(detectability), length(detectability) == 1,
    is.numeric(minSampleCoverage), length(minSampleCoverage) == 1
  )
  
  devMode <- getOption("HaplotypR.devel")
  if(is.null(devMode))
    devMode <- F
  
  if(is.null(referenceSequence)){
    referenceSequence <- markerTable$ReferenceSequence
    names(referenceSequence) <- markerTable$MarkerID
  }
    
  
  ## TODO

  # check if input file exists
  sampleTable <- sampleTable[!is.na(sampleTable$FileR1),]
  sampleTable <- sampleTable[file.exists(sampleTable$FileR1),]
  # filter read file
  newFileR1 <- file.path("filtered", basename(sampleTable$FileR1))
  dir.create("filtered")
  #file.remove(newFileR1)
  file.remove(list.files("tmp", full.names = T))
  numReads <- lapply(seq_along(sampleTable$FileR1), function(ii){
    sr <- readFastq(sampleTable$FileR1[ii])
    sr <- sr[0<width(sr)]
    #sr <- sr[0<width(sr) & width(sr)<2100]
    sr <- rmNreads(sr)
    writeFastq(sr, newFileR1[[ii]], compress = T)
  })
  sampleTable$FileR1 <- newFileR1
  sampleTable$numReads <- unlist(numReads)
  #sampleTable <- sampleTable[sampleTable$numReads>50,]
  #write.table(sampleTable, "demultiplexMarkerSummary_MinION_filt.txt", sep="\t", row.names=F)
  
  # run dada2
  err <- learnErrors(sampleTable$FileR1, multithread=F)
  # pdf("plotDataErrMinIon.pdf")
  # plotErrors(err, nominalQ=TRUE)
  # dev.off()
  dadas <- dada(sampleTable$FileR1, err=err, multithread=TRUE)
  seqtab <- makeSequenceTable(dadas)
  # remove chimera
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  #sum(seqtab.nochim)/sum(seqtab)
  seqtabLst <- split(as.data.frame(seqtab.nochim), sampleTable$MarkerID)
  #saveRDS(seqtabLst, file.path(".", "seqtabMinION.RDS"))
  
  selMarker <- unique(sampleTable$MarkerID)
  resultsLst <- lapply(unique(sampleTable$MarkerID), function(nm){
    # nm <- "MH_ama1_D2_18"
    message(nm)
    # get read counts and rename sample
    markTab <- seqtabLst[[nm]]
    markTab <- markTab[,colSums(markTab)>0, drop=F]
    rownames(markTab) <- sampleTable$SampleID[match(rownames(markTab), basename(sampleTable$FileR1))]
    # get haplotype sequence
    haplotyp <- DNAStringSet(colnames(markTab))
    names(haplotyp) <- paste(nm, "_", 1:length(haplotyp), sep="")
    colnames(markTab) <- names(haplotyp)
    
    # check sample cutoff
    totalReads <- rowSums(markTab, na.rm=T)
    # check minority clone cutoff
    minCov <- colSums(x, na.rm = T) * detectability
    minCov[minCov < minHaplotypCoverage] <- minHaplotypCoverage
    markTab[!minCov,] <- 0
    # check for length polymorphism resp indels
    refLength <- nchar(referenceSequence[nm])
    noIndel <- refLength==width(haplotyp)
    #noIndel <- (refLength-1)<width(haplotyp) & width(haplotyp)<(refLength+1) 
    markTab <- markTab[,noIndel, drop=F]

    
    
    # callHaplotype <- function(x, detectability=1/100, minHaplotypCoverage=3, minReplicate=2, minSampleCoverage=300)
    # haplotypesSample <- callHaplotype(haplotypesSample, minHaplotypCoverage=minHaplotypCoverage, minReplicate=minReplicate, 
    #                                   detectability=detectability, minSampleCoverage=1)
    


    keepHaplotype <- colSums(markTab)>0
    markTab <- markTab[,keepHaplotype, drop=F]
    

    # markTab <- cbind(markTab, Noise=rowSums(Noise1, na.rm=T)+rowSums(Noise2, na.rm=T))
    markTab <- cbind(markTab, Censored=totReads-rowSums(markTab, na.rm=T))
    Noise=totReads-rowSums(markTab, na.rm=T)
    rownames(markTab) <- selMH
    writeFasta(haplotyp, file=file.path(outputDir, sprintf("HaplotypeList_MinION_MH_CO_%s.fasta", nm)))
    


    #names(haplotyp) <- paste0(nm, "_MID", 1:length(haplotyp), sep="")
    #writeFasta(haplotyp, file=file.path(outputDir, sprintf("HaplotypeList_MinION_%s.fasta", nm)))
    #colnames(markTab) <- names(haplotyp)
    return(list(markTab, haplotyp))
  })
  names(resultsLst) <- selMarker

  
  #finalHaplotyopList <- writeHaplotypList(resultsLstMinMH)
  #write.csv(finalHaplotyopList, file=file.path(outputDir, "finalHaplotypList_MinION_MH.csv"), row.names=F)
  


  
  
  ## End TODO
  
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
                                                  minHaplotypCoverage, minSampleCoverage, minReplicate, detectability, marker, postfix)),
                sep="\t", row.names=F, col.names=T)
  
  
  # check replicates
  idx <- split(1:dim(haplotypesSample)[2], samTab$SampleID)
  markerRes <- lapply(idx, function(i){
    tab <- callHaplotype(haplotypesSample[,i, drop=F], minHaplotypCoverage=minHaplotypCoverage, 
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
}

