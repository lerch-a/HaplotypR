
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
  
  return(resTab[, c("SampleID","SampleName","BarcodePair","ReadNumbers","FileR1","FileR2")])
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
  res <- merge.data.frame(dePlexFiles[,c("SampleID", "SampleName","BarcodePair")], res, by="BarcodePair")
  
  write.table(res, file.path(out, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
  
  return(res)
}



runConcatReads <- function(input, output, session, volumes){

  #project <- createProject()
  project <- list()
  
  # Output dir
  out <- isolate(baseOutDir())
  if(is.null(out)){
    showNotification("Output Directory missing", closeButton = T, type="error")
    return(NULL)
  }
  project$outDir <- out
  
  markerTab <- isolate(primer())
  if(is.null(markerTab)){
    showNotification("Primer sequence missing", closeButton = T, type="error")
    return(NULL)
  }
  project$marker <- markerTab$MarkerID
  
  # Options
  wTrim <- isolate(input$withTrim)
  if(!wTrim){
    numNtF <- NULL
    numNtR <- NULL
    concatSuffix <- "_bind_full"
    project$concatSuffix <- concatSuffix
  }else{
    numNtF <- isolate(input$trim_F)
    numNtR <- isolate(input$trim_R)
    concatSuffix <- sprintf("_bind%.0f_%.0f", numNtF, numNtR)
    project$concatSuffix <- concatSuffix
    project$trimParameter <- c(Fwd=numNtF, Rev=numNtR)
    refSeq <- as.character(markerTab$ReferenceSequence)
    refSeq <- DNAStringSet(paste(substr(refSeq, 1,numNtF), substr(refSeq, nchar(refSeq)+1-numNtR, nchar(refSeq)), sep=""))
    names(refSeq) <- markerTab$MarkerID
    project$refSeq <- refSeq
    lapply(seq_along(refSeq), function(i){
      writeFasta(refSeq[i], file.path(out, paste(names(refSeq)[i], concatSuffix, ".fasta", sep="")))
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
  
  res <- cbind(dePlexFiles[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], res)
  write.table(res, file.path(out, sprintf("processedReadSummary%s.txt", concatSuffix)), sep="\t", row.names=F)
  
  #project <- createProject()
  project$readsFileLst <- res
  return(res)
}


runCallGenotype <- function(input, output, session, volumes){
  prj <- isolate(input$projectSelectedG)
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
  fn <- list.files(out, "processedReadSummary")
  fn <- fn[grep(prj, fn)]
  if(!file.exists(file.path(out, fn))){
    showNotification("Summary file of concatenate paired reads is missing. Did you run this step?", closeButton = T, type="error")
    return(NULL)
  }
  sampleTable <- read.delim(file.path(out, fn))
  #sampleTable <- isolate(outConcatReads())
  sampleTable <- sampleTable[sampleTable$MarkerID == marker,]
  sampleTable <- sampleTable[!(is.na(sampleTable$ReadFile) | is.na(sampleTable$SampleID)),]
  # Load Reference
  if(grepl("bind_full.fastq", sampleTable$ReadFile[1])){
    postfix <- "bind_full"
    markerTab <- isolate(primer())
    refSeq <- as.character(markerTab$ReferenceSequence)
    refSeq <- DNAStringSet(refSeq)
    names(refSeq) <- markerTab$MarkerID
    refSeq <- refSeq[grep(marker, markerTab$MarkerID)]
  }else{
    postfix <- sub("^.*_bind", "_bind", sampleTable$ReadFile[1])
    postfix <- sub("\\.fastq.*$", "", postfix)
    fn <- file.path(out, paste(marker, postfix, ".fasta", sep=""))
    if(!file.exists(fn)){
      showNotification("Reference sequence not found.", closeButton = T, type="error")
      return()
    }
    refSeq <- readDNAStringSet(fn)
  }
  
  # Configure shiny progress
  progress <- Progress$new(session, min=0, max=length(sampleTable$ReadFile))
  on.exit(progress$close())
  progress$set(message=sprintf('Calculation of mis-match rate in progress. This may take a while...'), detail=NULL, value=0)
  
  updateProgress <- function(detail=NULL, value=NULL){
    progress$set(detail=detail, value=value)
  }
  
  # Calculate mismatch rate
  seqErrLst <- calculateMismatchFrequencies(as.character(sampleTable$ReadFile), refSeq, minCoverage=100L, progressReport=updateProgress)
  names(seqErrLst) <- sampleTable$SampleID
  seqErr <- do.call(cbind, lapply(seqErrLst, function(l){
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  write.table(seqErr, file.path(out, sprintf("mismatchRate_rate%.0f_occ%i_%s%s.txt", 
                                             minMMrate*100, minOccGen, marker, postfix)), sep="\t", row.names=F)

  potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
  snpRef <- unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq, start=snp, width=1))
  }))
  snpLst <- cbind(Chr=names(refSeq), Pos=potSNP, Ref=snpRef, Alt="N")
  write.table(snpLst, file=file.path(out, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt", 
                                                  minMMrate*100, minOccGen, marker, postfix)), 
              row.names=F, col.names=T, sep="\t", quote=F)
  
  png(file.path(out, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png", 
                             minMMrate*100, minOccGen, marker, postfix)), 
      width=1500 , height=600)
  matplot(seqErr, type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1),
          ylab="Mismatch Rate", xlab="Base Position", main=marker, cex.axis=2, cex.lab=2)
  abline(v=snpLst[,"Pos"], lty=2, col="grey")
  abline(h=minMMrate, lty=1, col="red")
  dev.off()
#   return(list(seqErrors=seqErr, potentialSNP=potSNP))
# }
# 
# plotGenotype <- function(input, output, session, volumes){
  # genLst <- outCallGen()
  # if(is.null(genLst)) return(NULL)
  matplot(seqErr, type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1),
          ylab="Mismatch Rate", xlab="Base Position", cex.axis=1.5, cex.lab=1.5)
  abline(v=snpLst[,"Pos"], lty=2, col="grey")
  abline(h=minMMrate, lty=1, col="red")
}

runCreateHaplotypOverview <- function(input, output, session, volumes){
  prj <- isolate(input$projectSelectedH)
  marker <- isolate(input$markerSelectedH)
  if(is.null(marker)) return(NULL)
  if(marker=="none"){
    showNotification("Select Marker", closeButton = T, type="error")
    return(NULL)
  }

  # Output Dir
  out <- isolate(baseOutDir())
  
  # Load sample
  fn <- list.files(out, "processedReadSummary")
  fn <- fn[grep(prj, fn)]
  if(!file.exists(file.path(out, fn))){
    showNotification("Summary file of concatenate paired reads is missing. Did you run this step?", closeButton = T, type="error")
    return(NULL)
  }
  sampleTable <- read.delim(file.path(out, fn))
  #sampleTable <- isolate(outConcatReads())
  sampleTable <- sampleTable[sampleTable$MarkerID == marker,]
  sampleTable <- sampleTable[!is.na(sampleTable$ReadFile) & !is.na(sampleTable$SampleID),]
  sampleTable <- sampleTable[grep("NA", sampleTable$BarcodePair, invert=T),]
  if(grepl("_full.fastq", sampleTable$ReadFile[1])){
    postfix <- "full"
  }else{
    postfix <- sub("^.*_bind", "_bind", sampleTable$ReadFile[1])
    postfix <- sub("\\.fastq.*$", "", postfix)
  }
  
  # Output dir
  out <- isolate(baseOutDir())
  outFreqReads <- file.path(out, "frequencyFiles")
  if(!file.exists(outFreqReads)){
    dir.create(outFreqReads)
  }
  # SNP list
  if(input$fSNP){
	  potSNPLst <- outPotSNP()
	  if(is.null(potSNPLst)){
	    showNotification("Potential SNP list is missing. Did you run 'Call Genotype'?", closeButton = T, type="error")
	    return(NULL)
	  }
  }else{
  	potSNPLst <- NULL
  }
  
  # Configure shiny progress
  progressMain <- Progress$new(session, min=0, max=8)
  progressMain$set(message=sprintf('Call Haplotype is in progress. This may take a while...'), detail=NULL, value=1)
  progress <- Progress$new(session, min=0, max=length(sampleTable$ReadFile))
  on.exit({ 
    progress$close()
    progressMain$close()
  })
  progress$set(message=sprintf('Dereplication is in progress.'), detail=NULL, value=0)
  updateProgress <- function(detail=NULL, value=NULL){
    progress$set(detail=detail, value=value)
  }
  tab <- createContingencyTable(as.character(sampleTable$ReadFile), dereplicated=F, inputFormat="fastq", outputDir=outFreqReads,
                               sampleNames=as.character(sampleTable$SampleID), replicatNames=NULL, 
                               haplotypeUIDFile=NULL, progressReport=updateProgress)
  if(getAppOptionDevel())
    write.table(tab, file=file.path(out, sprintf("contingencyTable_%s%s.txt", marker, postfix)), sep="\t")
  
  fnAllSeq <- file.path(out, sprintf("allSequences_%s%s.fasta", marker, postfix))
  file.rename(file.path(outFreqReads, "allHaplotypes.fa"), fnAllSeq)

  frqfile <- file.path(outFreqReads, sub(".fastq.gz$", "_hapFreq.fa", basename(as.character(sampleTable$ReadFile))))
  #frqfile <- list.files(outFreqReads, pattern="hapFreq.fa$", full.names=T)
  prefix <- sub("_hapFreq\\.fa", "", basename(frqfile))
  
  outCluster <- file.path(out, "cluster")
  if(!file.exists(outCluster)){
    dir.create(outCluster)
  }
  outCluster <- file.path(outCluster, marker)
  if(!file.exists(outCluster)){
    dir.create(outCluster)
  }
  
  progressMain$set(value=2)
  progress$set(message=sprintf('Read clustering is in progress.'), detail=NULL, value=0)
  # run cluster
  clusterFilenames <- clusterReads(frqfile, outCluster, prefix, progressReport=updateProgress)

  # check chimeras
  progressMain$set(value=3)
  if(isolate(input$cChimera)){
    progress$set(message=sprintf('Check for chimera is in progress.'), detail=NULL, value=0)
    chimeraFilenames <- checkChimeras(clusterFilenames[,"RepresentativeFile"], method="vsearch", progressReport=updateProgress)
  }else{
    chimeraFilenames <- NULL
  }

  # check indels
  if(isolate(input$cIndels)){
    if(postfix=="full"){
      markerTab <- isolate(primer())
      refSeq <- as.character(markerTab$ReferenceSequence)
      refSeq <- DNAStringSet(refSeq)
      names(refSeq) <- markerTab$MarkerID
      refSeq <- refSeq[grep(marker, markerTab$MarkerID)]
    }else{
      refSeq <- sub("^.*_bind", "_bind", as.character(sampleTable$ReadFile)[1])
      refSeq <- sub("fastq.*$", "fasta", refSeq)
      refSeq <- file.path(out, paste(marker, refSeq, sep=""))
      refSeq <- readDNAStringSet(refSeq, format="fasta")
    }
  }else{
    refSeq <- NULL
  }

  progressMain$set(value=4)
  progress$set(message=sprintf('Create Haplotype Overview table in progress.'), detail="", value=length(sampleTable$ReadFile))
  overviewHap <- createHaplotypOverviewTable(fnAllSeq,
                                             clusterFilenames, chimeraFilenames,
                                             refSeq, potSNPLst)
  
  ## Final Haplotype
  progressMain$set(value=5)
  overviewHap$FinalHaplotype <- factor(NA, levels = c("Singelton", "Chimera", "Indels", marker))
  overviewHap[overviewHap$representatives, "FinalHaplotype"] <- marker
  overviewHap[overviewHap$singelton, "FinalHaplotype"] <- "Singelton"
  
  if(input$cIndels){
  	overviewHap[!is.na(overviewHap$indels) & overviewHap$indels & overviewHap$representatives, "FinalHaplotype"] <- "Indels"
  }
  
  if(input$cChimera){
  	#overviewHap[!is.na(overviewHap$chimera) & is.na(overviewHap$nonchimera), "FinalHaplotype"] <- "Chimera"
  	overviewHap[!is.na(overviewHap$chimeraScore) & is.na(overviewHap$nonchimeraScore), "FinalHaplotype"] <- "Chimera"
  }

  idx <- overviewHap$FinalHaplotype %in% marker
  if(input$fSNP){
	  snpLst <- unique(overviewHap$snps[idx])
	  names(snpLst) <- paste(marker, seq_along(snpLst), sep="-")
	  levels(overviewHap$FinalHaplotype) <- c(levels(overviewHap$FinalHaplotype), names(snpLst))
	  overviewHap$FinalHaplotype[idx] <- names(snpLst)[match(overviewHap$snps[idx], snpLst)]
  }else{
  	hapNames <- paste(marker, 1:sum(idx), sep="-")
  	levels(overviewHap$FinalHaplotype) <- c(levels(overviewHap$FinalHaplotype), hapNames)
  	overviewHap$FinalHaplotype[idx] <- hapNames
  }
  overviewHap <- as.data.frame(overviewHap)
  
  progressMain$set(value=6)
  write.table(overviewHap, file=file.path(out, sprintf("HaplotypeOverviewTable_%s%s.txt", marker, postfix)), sep="\t")
  
#   return(overviewHap)
# }
# 
# runCallHaplotype  <- function(input, output, session, volumes){

  # Options
  minCov <- isolate(input$minCoverage)
  maxSensitivity <- isolate(input$maxSensitivity)
  minOccHap <- isolate(input$minOccHap)
  minCovSample <- isolate(input$minCovSample)
  
  # marker <- isolate(input$markerSelectedH)
  # out <- isolate(baseOutDir())
  
  #####
  # Count table repfile
  progressMain$set(value=7)
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
  haplotypesSample <- haplotypesSample[rowSums(haplotypesSample)>0,]
  colnames(haplotypesSample) <- sub(".representatives.fasta", "", basename(repfile))

  overviewHap <- overviewHap[rownames(haplotypesSample),]

  # set final haplotype names
  rownames(haplotypesSample) <- overviewHap[rownames(haplotypesSample), "FinalHaplotype"]
  haplotypesSample <- reclusterHaplotypesTable(haplotypesSample)
  if(getAppOptionDevel())
    write.table(cbind(HaplotypNames=rownames(haplotypesSample),haplotypesSample), 
              file=file.path(out, sprintf("rawHaplotypeTable_%s%s.txt", marker, postfix)), sep="\t", row.names=F, col.names=T)
  #haplotypesSample <- reclusterHaplotypesTable(haplotypesSample)
  
  # Apply cut-off haplotype only in 1 sample
  
  haplotypesSample <- callHaplotype(haplotypesSample, minHaplotypCoverage=minCov, minReplicate=minOccHap, 
                                    detectability=maxSensitivity, minSampleCoverage=1, reportBackground=T)

  #haplotypesSample <- as.data.frame(cbind(HaplotypNames=rownames(haplotypesSample), haplotypesSample))
  #haplotypesSample <- reclusterHaplotypesTable(haplotypesSample)
  progressMain$set(value=7)
  if(getAppOptionDevel())
    write.table(cbind(HaplotypNames=rownames(haplotypesSample), haplotypesSample), 
  						file=file.path(out, sprintf("finalHaplotypeTable_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt", 
  						                            minCov, minCovSample, minOccHap, maxSensitivity, marker, postfix)), 
  													 sep="\t", row.names=F, col.names=T)

  # check replicates
  if(isolate(input$checkReplicates)){
    idx <- split(1:dim(haplotypesSample)[2], sampleTable$SampleName)
    lst <- lapply(idx, function(i){
      tab <- callHaplotype(haplotypesSample[,i, drop=F], minHaplotypCoverage=minCov, 
                           minReplicate=length(i), detectability=maxSensitivity, minSampleCoverage=minCovSample, 
                           reportBackground=T)
      tab <- cbind(sampleTable[rep(i,each=dim(tab)[1]), c("SampleID","SampleName","MarkerID")], 
                   Haplotype=rownames(tab), Reads=as.integer(tab))
      colnames(tab) <- c("SampleID","SampleName","MarkerID","Haplotype","Reads")
      rownames(tab) <- NULL
      #check individual samples for chimera
      tmpTab <- as.data.frame(tapply(tab$Reads, tab$Haplotype, sum))
      tmpTab$Haplotype <- rownames(tmpTab)
      colnames(tmpTab) <- c("Reads","Haplotype")
      hIdx <- grep(marker, tmpTab$Haplotype)
      chim <- NULL
      if(length(hIdx)>2){
        chim <- flagChimera(tmpTab[hIdx,], overviewHap)
      }
      rm(tmpTab)
      tab$FlagChimera <- tab$Haplotype %in% chim
      return(tab)
    })
  }else{
    lst <- lapply(1:dim(haplotypesSample)[2], function(i){
      tab <- callHaplotype(haplotypesSample[,i, drop=F], minHaplotypCoverage=minCov, 
                           minReplicate=1, detectability=maxSensitivity, minSampleCoverage=minCovSample, 
                           reportBackground=T)
      suppressWarnings(
        tab <- cbind(sampleTable[i,c("SampleID","SampleName","MarkerID")], 
                     Haplotype=rownames(tab), Reads=as.integer(tab))
      )
      colnames(tab) <- c("SampleID","SampleName","MarkerID","Haplotype","Reads")
      rownames(tab) <- NULL
      #check individual samples for chimera
      hIdx <- grep(marker, tab$Haplotype)
      chim <- NULL
      if(length(hIdx)>2){
        chim <- flagChimera(tab[hIdx,], overviewHap)
      }
      tab$FlagChimera <- tab$Haplotype %in% chim
      return(tab)
    })
  }
  lst <- do.call(rbind, lst)
  #TODO check individual samples for chimera

  write.table(lst, 
  						file=file.path(out, sprintf("finalHaplotypeList_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt", 
  						                            minCov, minCovSample, minOccHap, maxSensitivity, marker, postfix)), 
  						sep="\t", row.names=F, col.names=T)
  
  
  progressMain$set(value=8)
  #return(as.data.frame(cbind(HaplotypNames=rownames(haplotypesSample))))
  return(as.data.frame(lst))
}

flagChimera <- function(hapTable, overviewHap){
  snps <- overviewHap[overviewHap$FinalHaplotype %in% hapTable$Haplotype, c("snps", "FinalHaplotype")]
  snps <- snps[!duplicated(snps),]
  rownames(snps) <- snps$FinalHaplotype
  snps$snps <- as.character(snps$snps)
  snps <- snps[as.character(hapTable$Haplotype), "snps"]
  names(snps) <- hapTable$Haplotype
  snps <- snps[order(hapTable$Reads, decreasing=T)]
  chim <- findChimeras(snps)
  return(chim[,"chimera"])
}
