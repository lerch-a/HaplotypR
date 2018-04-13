
# compute coverage matrix C of size (N x M)
# rows correspond to N samples
# cols correspond to M haplotypes 
# c_ij is integer >= 0 for (i, j) in (N x M)  
createContingencyTable <- function(inputFiles, dereplicated=F, inputFormat="fasta", 
                                   outputDir=".", sampleNames, replicatNames=NULL, freqSplitPattern="_", 
                                   haplotypeUIDFile=NULL, progressReport=message, include_seq=FALSE, verbose=FALSE) {
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
    if (!dereplicated) {
      readFreq <- tables(sread(inputReads), n=NULL)$top
      haplotypes <- DNAStringSet(names(readFreq))
      names(readFreq) <- NULL
    } else {
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
  #print(env$allHaplotypes)
  writeFasta(env$allHaplotypes, file.path(outputDir, "allHaplotypes.fa"))
  
  cat("\nbuilding contigency table...\n")
  freq <- integer(length(env$allHaplotypes))
  names(freq) <- names(env$allHaplotypes)
  contingencyTable <- lapply(contingencyList, function(l) {
    #idx <- match(l$UID, names(env$allHaplotypes))
    freq[as.character(l$UID)] <- l$Frequency
    return(freq)
  })
  contingencyTable <- do.call(cbind, contingencyTable)
  #print(as.character(env$allHaplotypes))
  rownames(contingencyTable) <- names(env$allHaplotypes)
  colnames(contingencyTable) <- names(contingencyList) 
  contingencyTable = as.data.frame(contingencyTable)
  if (include_seq) 
    contingencyTable[["seq"]] <- as.character(env$allHaplotypes)
  
  return(contingencyTable)
}

# call haplotype, if passes all filters 
callHaplotype <- function(x, detectability=1/100, minHaplotypCoverage=3, minReplicate=2, 
                          minSampleCoverage=300, reportBackground=T, defineBackground=NULL, ...) {

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

# compute final table of called haplotypes 
createFinalHaplotypeTable <- function(outputDir, sampleTable, markerTab, snpLst, refSeq, postfix, 
                                      minHaplotypCoverage=2, minReplicate=3, detectability=0.01, minSampleCoverage=25, 
                                      include_seq=FALSE, verbose=FALSE) {
    source("/Users/tfarrell/Tools/HaplotypR/R/clusterFunction.R")
    source("/Users/tfarrell/Tools/HaplotypR/R/filterHaplotype.R")
  outFreqFiles <- file.path(outputDir, "haplotype_freq_files")
  dir.create(outFreqFiles)
  res <- lapply(markerTab$MarkerID, function(marker){
    samTab <- sampleTable[sampleTable$MarkerID == marker,]
    potSNPLst <- snpLst[[marker]]
    prefix <- sub(".fastq.gz$", "", basename(as.character(samTab$ReadFile)))
    
    # Create frequency files and count table
    contingency_file = file.path(outputDir, sprintf("contingencyTable_%s%s.txt", marker, postfix))
    if (!file.exists(contingency_file)) { 
        if (verbose) cat("\ncomputing contingency table...\n")
        tab <- createContingencyTable(as.character(samTab$ReadFile), dereplicated=F, 
                                      inputFormat="fastq", outputDir=outFreqFiles,
                                      sampleNames=as.character(samTab$SampleID), replicatNames=NULL, 
                                      haplotypeUIDFile=NULL, include_seq=include_seq, verbose=verbose)
        if (verbose) cat("\nwriting contingency table...\n")
        write.table(tab, file=contingency_file, sep="\t")
    } else { 
        if (verbose) cat("\nreading in contingency table...\n")
        #tab = read.table(contingency_file) 
    }
    fnAllSeq <- file.path(outFreqFiles, sprintf("allSequences_%s%s.fasta", marker, postfix))
    file.rename(file.path(outFreqFiles, "allHaplotypes.fa"), fnAllSeq)
    frqfile <- file.path(outFreqFiles, paste(prefix, "_hapFreq.fa", sep=""))
    
    # run cluster with Rswarm package
    if (verbose) cat("\nrunning Rswarm...\n")
    outCluster <- file.path(outputDir, "cluster", marker)
    if (!dir.exists(outCluster)) dir.create(outCluster, recursive=T)
    clusterFilenames <- clusterReads(frqfile, outCluster, prefix)
    
    # check for chimeras with Rvsearch package
    if (verbose) cat("\nrunning Rsearch...\n")
    chimeraFilenames <- checkChimeras(clusterFilenames[,"RepresentativeFile"], method="vsearch")
    
    # create an overview table
    if (verbose) cat("\ncomputing final haplotype table...\n")
    overviewHap <- createHaplotypOverviewTable(fnAllSeq,
                                               clusterFilenames, chimeraFilenames,
                                               refSeq[marker], potSNPLst, verbose=verbose)
    ## create final haplotype
    overviewHap$FinalHaplotype <- factor(NA, levels = c("Singelton", "Chimera", "Indels", as.character(marker)))
    overviewHap[overviewHap$representatives, "FinalHaplotype"] <- as.character(marker)
    # Set singelton
    overviewHap[overviewHap$singelton, "FinalHaplotype"] <- "Singelton"
    # Set Indel
    overviewHap[!is.na(overviewHap$indels) & overviewHap$indels & overviewHap$representatives, "FinalHaplotype"] <- "Indels"
    # Set chimera
    overviewHap[!is.na(overviewHap$chimeraScore) & is.na(overviewHap$nonchimeraScore), "FinalHaplotype"] <- "Chimera"
    #print(overviewHap)
    # Cluster identical SNPs pattern
    idx <- overviewHap$FinalHaplotype %in% marker
    #print(idx)
    snps <- unique(overviewHap$snps[idx])
    #print(snps)
    names(snps) <- paste(marker, seq_along(snps), sep="-")
    levels(overviewHap$FinalHaplotype) <- c(levels(overviewHap$FinalHaplotype), names(snps))
    overviewHap$FinalHaplotype[idx] <- names(snps)[match(overviewHap$snps[idx], snps)]
    overviewHap <- as.data.frame(overviewHap)
    #write.table(overviewHap, file=file.path(outputDir, sprintf("HaplotypeOverviewTable_%s%s.txt", marker, postfix)), sep="\t")
    
    # runCallHaplotype  <- function(input, output, session, volumes){
    
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
    #write.table(cbind(HaplotypNames=rownames(haplotypesSample),haplotypesSample), 
    #              file=file.path(outputDir, sprintf("rawHaplotypeTable_%s%s.txt", marker, postfix)), sep="\t", row.names=F, col.names=T)
    haplotypesSample <- reclusterHaplotypesTable(haplotypesSample)
    #write.table(cbind(HaplotypNames=rownames(haplotypesSample),haplotypesSample), 
    #              file=file.path(outputDir, sprintf("reclusteredHaplotypeTable_%s%s.txt", marker, postfix)), sep="\t", row.names=F, col.names=T)
    
    
    # Apply cut-off haplotype only in 1 sample
    haplotypesSample <- callHaplotype(haplotypesSample, minHaplotypCoverage=minHaplotypCoverage, minReplicate=minReplicate, 
                                      detectability=detectability, minSampleCoverage=1)
    
    # write.table(cbind(HaplotypNames=rownames(haplotypesSample), haplotypesSample), 
    #             file=file.path(outputDir, sprintf("finalHaplotypeTable_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt", 
    #                                         minCov, minCovSample, minOccHap, detectionLimit, marker, postfix)), 
    #             sep="\t", row.names=F, col.names=T)
    
    # check replicates
    idx <- split(1:dim(haplotypesSample)[2], samTab$SampleName)
    lst <- lapply(idx, function(i){
      tab <- callHaplotype(haplotypesSample[,i, drop=F], minHaplotypCoverage=minHaplotypCoverage, 
                           minReplicate=length(i), detectability=detectability, minSampleCoverage=minSampleCoverage, 
                           reportBackground=T)
      tab <- cbind(samTab[rep(i,each=dim(tab)[1]), c("SampleID","SampleName","MarkerID")], 
                   Haplotype=rownames(tab), Reads=as.integer(tab))
      colnames(tab) <- c("SampleID","SampleName","MarkerID","Haplotype","Reads")
      rownames(tab) <- NULL
      #check individual samples for chimera
      if(length(i)){
        tmpTab <- as.data.frame(tapply(tab$Reads, tab$Haplotype, sum))
        tmpTab$Haplotype <- rownames(tmpTab)
        colnames(tmpTab) <- c("Reads","Haplotype")
        hIdx <- grep(marker, tmpTab$Haplotype)
      }else{
        hIdx <- grep(marker, tab$Haplotype)
      }
      chim <- NULL
      if(length(hIdx)>2){
        chim <- flagChimera(tmpTab[hIdx,], overviewHap)
      }
      rm(tmpTab)
      #tab$FlagChimera <- tab$Haplotype %in% chim
      return(tab)
    })
    lst <- do.call(rbind, lst)
    
    
    write.table(lst, 
                file=file.path(outputDir, sprintf("finalHaplotypeList_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt", 
                                                  minHaplotypCoverage, minSampleCoverage, minReplicate, detectability, 
                                                  as.character(marker), postfix)), 
                sep="\t", row.names=F, col.names=T)
    rownames(lst) <- NULL
    return(lst)
  })
  names(res) <- markerTab$MarkerID
  return(res)
}
