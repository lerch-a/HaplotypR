
createSampleFile <- function(input, output, session, volumes){
  fnB1 <- as.character(isolate(parseFilePaths(volumes, input$fileB1)$datapath))
  fnB2 <- as.character(isolate(parseFilePaths(volumes, input$fileB2)$datapath))
  if(is.null(fnB1)){
    showNotification("Forwarde barcode file missing.", closeButton = T, type="error")
    return(NULL)
  }
  if(is.null(fnB2)){
    showNotification("Reverse barcode file missing.", closeButton = T, type="error")
    return(NULL)
  }
  barcodeF <- readFasta(fnB1)
  barcodeR <- readFasta(fnB2)
  lst <- do.call(rbind, lapply(as.character(id(barcodeR)), function(br){
    data.frame(BarcodeID_F=as.character(id(barcodeF)), BarcodeID_R=br)
  }))
  lst <- as.data.frame(cbind(SampleID="", lst, AddInfo1="", AddInfo2=""))
  return(lst)
}

runDemultiplex <- function(input, output, session, volumes){
  
  fnR1 <- sub("/\\./", "/", file.path(volumes, isolate(fileR1())))
  fnR2 <- sub("/\\./", "/", file.path(volumes, isolate(fileR2())))
  fnB1 <- sub("/\\./", "/", file.path(volumes, isolate(fileB1())))
  fnB2 <- sub("/\\./", "/", file.path(volumes, isolate(fileB2())))
  wSampleName <- isolate(input$withSampleName)

  # read and check input values
  if(is.null(fnR1)){
    showNotification("Forward seuquence file missing.", closeButton = T, type="error")
    return(NULL)
  }
  if(is.null(fnR2)){
    showNotification("Reverse seuquence file missing.", closeButton = T, type="error")
    return(NULL)
  }
  if(is.null(fnB1)){
    showNotification("Forwarde barcode file missing.", closeButton = T, type="error")
    return(NULL)
  }
  if(is.null(fnB2)){
    showNotification("Reverse barcode file missing.", closeButton = T, type="error")
    return(NULL)
  }
  if(wSampleName){
    sampleTab <- isolate(sample())
    if(length(sampleTab)==0){
      showNotification("Sample file missing.", closeButton = T, type="error")
      return(NULL)
    }
  }
  
  out <- isolate(baseOutDir())
  
  # Demultiplexed sample files
  outDir <- file.path(out, "demultiplexBySample")
  if(!file.exists(outDir)){
    dir.create(outDir)
  }else{
    showNotification("Output dir exists", closeButton = T, type="error")
    return(NULL)
  }

  # Configure shiny progress
  progress <- Progress$new(session, min=0, max=1)
  on.exit(progress$close())
  progress$set(message='De-Multiplexing in progress. This may take a while...', value=1)
  updateProgress <- function(detail=NULL){
    progress$set(detail=detail)
  }
  
  # run demultiplexing
  res <- demultiplexReads(fnR1, fnR2, fnB1, fnB2, outDir, progressReport=updateProgress)

  if(wSampleName && length(sampleTab)!=0){
    res <- renameDemultiplexedFiles(sampleTab, res)
  }else{
    res <- cbind(SampleID=NA, res)
  }
  
  write.table(res, file.path(out, "demultiplexSampleSummary.txt"), sep="\t", row.names=F)
  
  return(res) 
}


renameDemultiplexedFiles <- function(sampleTab, resTab){

  sampleTab$BarcodePair <- paste(sampleTab$BarcodeID_F, sampleTab$BarcodeID_R, sep="-")
  resTab <- merge.data.frame(sampleTab, resTab, by="BarcodePair", all.y=T)

  resTab <- lapply(1:dim(resTab)[1], function(i){
    #if(!is.na(resTab$SampleID[i])){
      old <- as.character(resTab[i, c("FileR1","FileR2")])
      #new <- sub(resTab$BarcodePair[i], resTab$SampleID[i], old)
      new <- file.path(dirname(old), paste(resTab$SampleID[i], "_BC_", basename(old), sep=""))
      resTab[i, c("FileR1","FileR2")] <- new
      file.rename(as.character(old), as.character(new))
    #}
    return(resTab[i,])
  })
  resTab <- do.call(rbind, resTab)
  
  return(resTab[, c("SampleID","BarcodePair","ReadNumbers","FileR1","FileR2")])
}


runTruncatePrimer <- function(input, output, session, volumes){
  
  # Options
  numMM <- isolate(input$numMisMatch)
  wIndel <- isolate(input$withIndel)
  
  # Primer sequence
  markerTab <- isolate(primer())
  if(is.null(markerTab)){
    showNotification("Primer sequence missing", closeButton = T, type="error")
    return(NULL)
  }
  
  # Demultiplexed sample files
  out <- isolate(baseOutDir())
  dePlexFiles <- isolate(outDePlexSample())
  
  # Output dir
  outDir <- file.path(out, "demultiplexByMarker")
  if(!file.exists(outDir)){
    dir.create(outDir)
  }else{
    # showModal(modalDialog(
    #     helpText(sprintf("Output directory '%s' exists! Press OK to overwrite.", outDir)),
    #     title="Warning!",
    #     size="l",
    #     footer = tagList(
    #       modalButton("Cancel"),
    #       actionButton("ok", "OK")
    #     )
    #   ))
    showNotification("Output dir exists", closeButton = T, type="error")
    return(NULL)
  }

  # Configure shiny progress
  progress <- Progress$new(session, min=0, max=length(dePlexFiles$FileR1))
  on.exit(progress$close())
  updateProgress <- function(detail=NULL){
    progress$set(detail=detail)
  }
  
  # Run
  res <- lapply(1:dim(markerTab)[1], function(j){
    mID <- as.character(markerTab[j, "MarkerID"])
    adapterF <- as.character(markerTab[j, "Forward"])
    adapterR <- as.character(markerTab[j, "Reverse"])
    #dir.create(file.path(outDir, mID), showWarnings=F)

    res <- lapply(seq_along(dePlexFiles$FileR1), function(i){
      progress$set(message=sprintf('De-Multiplexing of marker %s in progress. This may take a while...', mID), detail=NULL, value=i)
      outputFile <- file.path(outDir, sub("R1\\.fastq.gz", mID, basename(as.character(dePlexFiles$FileR1)[i])))
      removePrimer(as.character(dePlexFiles$FileR1)[i], as.character(dePlexFiles$FileR2)[i], outputFile, 
                   adapterF, adapterR, max.mismatch=numMM, with.indels=wIndel, progressReport=updateProgress)
    })
    cbind(BarcodePair=as.character(dePlexFiles$BarcodePair), MarkerID=mID, do.call(rbind, res))
  })
  res <- do.call(rbind, res)
  res <- merge.data.frame(dePlexFiles[,c("SampleID", "BarcodePair")], res, by="BarcodePair")
  
  write.table(res, file.path(out, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
  
  return(res)
}



runConcatReads <- function(input, output, session, volumes){

  # Output dir
  out <- isolate(baseOutDir())
  if(is.null(out)){
    showNotification("Output Directory missing", closeButton = T, type="error")
    return(NULL)
  }
  
  # Options
  wTrim <- isolate(input$withTrim)
  if(!wTrim){
    numNtF <- NULL
    numNtR <- NULL
  }else{
    numNtF <- isolate(input$trim_F)
    numNtR <- isolate(input$trim_R)
    markerTab <- isolate(primer())
    if(is.null(markerTab)){
      showNotification("Primer sequence missing", closeButton = T, type="error")
      return(NULL)
    }
    refSeq <- as.character(markerTab$ReferenceSequence)
    refSeq <- DNAStringSet(paste(substr(refSeq, 1,numNtF), substr(refSeq, nchar(refSeq)+1-numNtR, nchar(refSeq)), sep=""))
    names(refSeq) <- markerTab$MarkerID
    lapply(seq_along(refSeq), function(i){
      writeFasta(refSeq[i], file.path(out, paste(names(refSeq)[i], "_bind", input$trim_F, "_", input$trim_R, ".fasta", sep="")))
    })
  }
  
  # Demultiplexed sample files
  outDePlexMarker <- file.path(out, "demultiplexMarkerSummary.txt")
  #outDePlexMarker <- outDePlexMarker()
  if(file.exists(outDePlexMarker)){
    dePlexFiles <- read.delim(outDePlexMarker)
  }else{
    showNotification("Input missing", closeButton = T, type="error")
    return(NULL)
  }
  dePlexFiles <- dePlexFiles[!is.na(dePlexFiles$FileR1),]
  
  # Output dir
  outDir <- file.path(out, "processedReads")
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  # else{
  #   showNotification("Output dir exists", closeButton = T, type="error")
  #   return(NULL)
  # }
  
  # Configure shiny progress
  progress <- Progress$new(session, min=0, max=length(dePlexFiles$FileR1))
  on.exit(progress$close())
  progress$set(message=sprintf('Concatenation of reads in progress. This may take a while...'), detail=NULL, value=0)
  
  updateProgress <- function(detail=NULL, value=NULL){
    progress$set(detail=detail, value=value)
  }
  
  # Run
  #res <- lapply(seq_along(dePlexFiles$FileR1), function(i){
    #progress$set(message=sprintf('Concatenation of reads in progress. This may take a while...'), detail=NULL, value=i)
    #outputFile <- file.path(outDir, mID, sub("R1\\.fastq.gz", mID, basename(as.character(dePlexFiles$FileR1)[i])))
  res <- bindAmpliconReads(as.character(dePlexFiles$FileR1), as.character(dePlexFiles$FileR2), outDir, 
                      read1Length=numNtF, read2Length=numNtR, progressReport=updateProgress)
  #})
#    cbind(BarcodePair=as.character(dePlexFiles$BarcodePair), MarkerID=mID, do.call(rbind, res))
  #res <- do.call(rbind, res)
  
  res <- cbind(dePlexFiles[,c("SampleID", "BarcodePair", "MarkerID")], res)
  
  write.table(res, file.path(out, "processedReadSummary.txt"), sep="\t", row.names=F)
  
  return(res)
}


runCallGenotype <- function(input, output, session, volumes){
  marker <- isolate(input$markerSelectedG)
  if(is.null(marker)) return(NULL)
  if(marker=="none"){
    showNotification("Select Marker", closeButton = T, type="error")
    return(NULL)
  }
  # Output Dir
  out <- isolate(baseOutDir())
  # Options
  minMMrate <- isolate(input$minMMrate)
  minOccGen <- isolate(input$minOccGen)
  # Load Sample   
  sampleTable <- isolate(outConcatReads())
  sampleTable <- sampleTable[sampleTable$MarkerID == marker,]
  sampleTable <- sampleTable[!(is.na(sampleTable$ReadFile) | is.na(sampleTable$SampleID)),]
  # Load Reference
  postfix <- sub("^.*_bind", "_bind", sampleTable$ReadFile[1])
  postfix <- sub("\\.fastq.*$", ".fasta", postfix)
  fn <- file.path(out, paste(marker, postfix, sep=""))
  if(!file.exists(fn)){
    showNotification("Reference sequence not found.", closeButton = T, type="error")
    return()
  }
  refSeq <- readDNAStringSet(fn)
  
  # Configure shiny progress
  progress <- Progress$new(session, min=0, max=length(sampleTable$ReadFile))
  on.exit(progress$close())
  progress$set(message=sprintf('Calculation of mis-match rate in progress. This may take a while...'), detail=NULL, value=0)
  
  updateProgress <- function(detail=NULL, value=NULL){
    progress$set(detail=detail, value=value)
  }
  
  # Calculate mismatch rate
  seqErrLst <- calculateMismatchFrequencies(as.character(sampleTable$ReadFile), refSeq, minCoverage=50L, progressReport=updateProgress)
  names(seqErrLst) <- sampleTable$SampleID
  seqErr <- do.call(cbind, lapply(seqErrLst, function(l){
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  write.table(seqErr, file.path(out, sprintf("mismatchRate_%s.txt", marker)), sep="\t", row.names=F)
  
  potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
  snpRef <- unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq, start=snp, width=1))
  }))
  snpLst <- cbind(Chr=names(refSeq), Pos=potSNP, Ref=snpRef, Alt="N")
  write.table(snpLst, file=file.path(out, sprintf("potentialSNPlist_rate%.0f_occ%i_%s.txt", 
                                                  minMMrate*100, minOccGen, marker)), 
              row.names=F, col.names=T, sep="\t", quote=F)
  
  return(list(seqErrors=seqErr, potentialSNP=potSNP))
}

plotGenotype <- function(input, output, session, volumes){
  genLst <- outCallGen()
  if(is.null(genLst)) return(NULL)
  matplot(genLst$seqErrors, type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1),
          ylab="Mismatch Rate", xlab="Base Position", cex.axis=1.5, cex.lab=1.5)
  abline(v=genLst$potentialSNP$Pos, lty=2, col="grey")
}

runCallHaplotype <- function(input, output, session, volumes){
  browser()
  marker <- isolate(input$markerSelectedH)
  if(is.null(marker)) return(NULL)
  if(marker=="none"){
    showNotification("Select Marker", closeButton = T, type="error")
    return(NULL)
  }
  # Options
  minMMrate <- isolate(input$minCoverage)
  maxSensitivity <- isolate(input$maxSensitivity)
  minOccGen <- isolate(input$minOccHap)
  # Load sample
  sampleTable <- isolate(outConcatReads())
  sampleTable <- sampleTable[sampleTable$MarkerID == marker,]
  sampleTable <- sampleTable[!(is.na(sampleTable$ReadFile) | is.na(sampleTable$SampleID)),]
  # Output dir
  out <- isolate(baseOutDir())
  outFreqReads <- file.path(out, "frequencyFiles")
  if(!file.exists(outFreqReads)){
    dir.create(outFreqReads)
  }
  # SNP list
  potSNPLst <- outPotSNP()
  if(is.null(potSNPLst)){
    showNotification("Potential SNP list is missing", closeButton = T, type="error")
    return(NULL)
  }

  tab <- createContingencyTable(as.character(sampleTable$FileR), dereplicated=F, inputFormat="fastq", outputDir=outFreqReads,
                               sampleNames=as.character(sampleTable$SampleID), replicatNames="", 
                               haplotypeUIDFile=NULL)
  browser()
  write.table(tab, file=file.path(out, "contingencyTable.txt"), sep="\t")
  file.rename(file.path(outFreqReads, "allHaplotypes.fa"), file.path(out, "allSequences.fasta"))

  frqfile <- list.files(outFreqReads, pattern="hapFreq.fa$", full.names=T)
  prefix <- sub("_hapFreq\\.fa", "", basename(frqfile))
  
  outCluster <- file.path(out, "cluster")
  if(!file.exists(outCluster)){
    dir.create(outCluster)
  }
  
  # run cluster
  clusterFilenames <- clusterReads(frqfile, outCluster, prefix)
  
  # check chimeras
  if(isolate(input$cChimera)){
    chimeraFilenames <- checkChimeras(clusterFilenames[,"RepresentativeFile"], method="vsearch")
  }else{
    chimeraFilenames <- NULL
  }

  # check indels
  if(isolate(input$cIndels)){
    refSeq <- sub("^.*_bind", "_bind", as.character(sampleTable$FileR)[1])
    refSeq <- sub("fastq.*$", "fasta", refSeq)
    refSeq <- file.path(out, paste(marker, refSeq, sep="_"))
    refSeq <- readDNAStringSet(refSeq, format="fasta")
  }else{
    refSeq <- NULL
  }

  overviewHap <- createHaplotypOverviewTable(file.path(out, "allSequences.fasta"),
                                             clusterFilenames, chimeraFilenames,
                                             refSeq, potSNPLst)
  
  ## Final Haplotype
  overviewHap$FinalHaplotype <- factor(NA, levels = c("Singelton", "Chimera", "Indels", marker))
  overviewHap[overviewHap$representatives, "FinalHaplotype"] <- marker
  overviewHap[overviewHap$singelton, "FinalHaplotype"] <- "Singelton"
  overviewHap[!is.na(overviewHap$chimera) & is.na(overviewHap$nonchimera), "FinalHaplotype"] <- "Chimera"
  overviewHap[!is.na(overviewHap$indels) & overviewHap$indels & overviewHap$representatives, "FinalHaplotype"] <- "Indels"

  
  
  
  # idx <- overviewHap$FinalHaplotype %in% marker
  # snpLst <- unique(overviewHap$snps[idx])
  # names(snpLst) <- paste(marker, seq_along(snpLst), sep="-")
  # levels(overviewHap$FinalHaplotype) <- c(levels(overviewHap$FinalHaplotype), names(snpLst))
  # overviewHap$FinalHaplotype[idx] <- names(snpLst)[match(overviewHap$snps[idx], snpLst)]
  write.table(as.data.frame(overviewHap), file=file.path(out, "HaplotypeOverviewTable.txt"), sep="\t")
  
  # #write.table(overviewHap, file=file.path(outDir, "finalHaplotypeOverviewTable.txt"), sep="\t")
  # overviewHap <- read.delim(file=file.path(outDir, "finalHaplotypeOverviewTable.txt"))
  # 
  # 
  # #####
  # # Count table repfile
  # 
  # repfile <- list.files(file.path(outDir, "cluster"), pattern="^[RDX].*representatives.fasta", full.names=T)
  # 
  # overviewHap <- read.delim(file.path(outDir, "finalHaplotypeOverviewTable.txt"))
  # 
  # hab <- rep(0, dim(overviewHap)[1])
  # names(hab) <- rownames(overviewHap)
  # 
  # haplotypesSample <- lapply(seq_along(repfile), function(i){
  #   sr1 <- readFasta(repfile[i])
  #   vec <- do.call(rbind, strsplit(as.character(id(sr1)), "_"))
  #   clusterResp <- vec[,1]
  #   clusterSize <- as.integer(vec[,2])
  #   hab[clusterResp] <- clusterSize
  #   hab
  # })
  # haplotypesSample <- do.call(cbind, haplotypesSample)
  # haplotypesSample <- haplotypesSample[rowSums(haplotypesSample)>0,]
  # colnames(haplotypesSample) <- sub(".representatives.fasta", "", basename(repfile))
  # dim(haplotypesSample)
  # 
  # overviewHap <- overviewHap[rownames(haplotypesSample),]
  # dim(overviewHap)
  # head(overviewHap)
  # table(overviewHap$representatives)
  # table(overviewHap[, "FinalHaplotype"])
  # 
  # # set final haplotype names
  # rownames(haplotypesSample) <- overviewHap[rownames(haplotypesSample), "FinalHaplotype"]
  # 
  # # recluster
  # dim(haplotypesSample)
  # haplotypesSample <- reclusterHaplotypesTable(haplotypesSample)
  # dim(haplotypesSample)
  # 

  data.frame(Test=1)
}
