

checkChimeras <- function(representativesFile, method="vsearch", progressReport=message) {
    # sort representatives 
    sortfile <- unlist(lapply(seq_along(representativesFile), function(i) {
        sortfile <- sub(".fasta", "_sort.fasta", representativesFile[i])
        if (!file.exists(sortfile)) { 
            #vsearch --sortbysize R-1_Rep1.representatives.fasta --output R-1_Rep1.representatives_sortbysize.fasta
            sr1 <- readFasta(representativesFile[i])
            clusterSize <- as.integer(do.call(rbind, strsplit(as.character(id(sr1)), "_"))[,2])
            sr1 <- sr1[order(clusterSize, decreasing=T)]
            writeFasta(ShortRead(sread(sr1), BStringSet(sub("_", ";size=", as.character(id(sr1))))), sortfile)
        } 
        return(sortfile)
    }))
    
    nonchimerafile <- sub("_sort.fasta", "_nonchimera.fasta", sortfile)
    chimerafile <- sub("_sort.fasta", "_chimera.fasta", sortfile)
    borderfile <- sub("_sort.fasta", "_bordchimera.fasta", sortfile)
    resfile <- sub("_sort.fasta", "_chimeraResults.txt", sortfile)
    # vsearch --uchime_denovo R-1_Rep1.representatives_nonchimeras.fasta --nonchimeras R-1_Rep1.representatives_nonchimeras2.fasta
    syscall <- paste("--uchime_denovo ", sortfile,"--mindiffs", 3, "--minh", 0.2, "--nonchimeras", nonchimerafile, 
                     "--chimeras", chimerafile, "--borderline", borderfile, "--uchimeout", resfile, sep=" ")
    
    # check and set progress report function
    if(!is.function(progressReport))
        progressReport <- message
    
    require(Rvsearch)
    files = c(nonchimerafile, chimerafile, borderfile, resfile)
    #Rvsearch:::.vsearchBin(args=syscall)
    lapply(seq_along(syscall), function(i) { 
        if (!all(sapply(files, file.exists))) {
            msg <- paste("Processing file", basename(representativesFile[i]), "...", sep=" ")
            progressReport(detail=msg, value=i)
            Rvsearch:::.vsearchBin(args=syscall[i])
        } else { 
            msg <- paste0("Retrieving existing files for ", basename(representativesFile[i]), "...")
            progressReport(detail=msg, value=i)
        }
    })
    
    return(cbind(NonchimeraFile=nonchimerafile, ChimeraFile=chimerafile, BorderchimeraFile=borderfile, ChimeraResultsFile=resfile))
}


createHaplotypOverviewTable <- function(allHaplotypesFilenames, clusterFilenames, chimeraFilenames, referenceSequence=NULL, 
                                        snpList, verbose=FALSE){
  
  sr1 <- readFasta(allHaplotypesFilenames)
  overviewHap <- data.frame(HaplotypesName=as.character(id(sr1)))
  rownames(overviewHap) <- overviewHap$HaplotypesName

  ######
  ## chimera type
  resfile <- chimeraFilenames[,"ChimeraResultsFile"]
  #sampleName <- sub(".representatives_chimeraResults.fasta", "", basename(resfile))
  haplotypes <- lapply(seq_along(resfile), function(i) {
    if (verbose) cat(paste0("Retrieving ", resfile[i], "...\n"))
    tryCatch({
            res <- read.delim(resfile[i], header = F)
            chimScore <- res[,1]
            vec <- do.call(rbind, strsplit(as.character(res[,2]), ";size="))
            clusterResp <- vec[,1]
            clusterSize <- as.integer(vec[,2])
            clusterResp <- clusterResp[clusterSize>1]
            chimScore <- chimScore[clusterSize>1]
            return(list(Chimera=clusterResp[chimScore>0], NonChimera=clusterResp[chimScore==0])) }, 
        warning=function (w) {  }, 
        error=function(e) { return(list(Chimera=character(0), NonChimera=character(0))) })
  })
  #names(haplotypes) <- sampleName
  chim <- table(unlist(lapply(haplotypes, "[[", 1)))
  nochim <- table(unlist(lapply(haplotypes, "[[", 2)))
  
  overviewHap[names(chim),"chimeraScore"] <- chim
  overviewHap[names(nochim),"nonchimeraScore"] <- nochim
  
  chimerafile <- chimeraFilenames[,"ChimeraFile"]
  sampleName <- sub(".representatives_chimera.fasta", "", basename(chimerafile))
  #sampleName <- do.call(rbind, strsplit(basename(chimerafile), "\\."))[,1]
  haplotypes <- lapply(seq_along(chimerafile), function(i) {
    sr1 <- readFasta(chimerafile[i])
    vec <- do.call(rbind, strsplit(as.character(id(sr1)), ";size="))
    clusterResp <- vec[,1]
    clusterSize <- as.integer(vec[,2])
    clusterResp <- clusterResp[clusterSize>1]
    if(length(clusterResp) < 0)
      clusterResp <- character(0)
    return(clusterResp)
  })
  names(haplotypes) <- sampleName
  tab <- table(unlist(haplotypes))
  overviewHap[names(tab),"chimera"] <- tab

  bordchimerafile <- chimeraFilenames[,"BorderchimeraFile"]
  haplotypes <- lapply(seq_along(bordchimerafile), function(i){
    sr1 <- readFasta(bordchimerafile[i])
    vec <- do.call(rbind, strsplit(as.character(id(sr1)), ";size="))
    clusterResp <- vec[,1]
    clusterSize <- as.integer(vec[,2])
    clusterResp <- clusterResp[clusterSize>1]
    if(length(clusterResp) < 0)
      clusterResp <- character(0)
    return(clusterResp)
  })
  names(haplotypes) <- sampleName
  tab <- table(unlist(haplotypes))
  overviewHap[names(tab),"bordchimera"] <- tab
  
  nonchimerafile <- chimeraFilenames[,"NonchimeraFile"]
  haplotypes <- lapply(seq_along(nonchimerafile), function(i) {
    sr1 <- readFasta(nonchimerafile[i])
    vec <- do.call(rbind, strsplit(as.character(id(sr1)), ";size="))
    clusterResp <- vec[,1]
    clusterSize <- as.integer(vec[,2])
    clusterResp <- clusterResp[clusterSize>1]
    if(length(clusterResp) < 0)
      clusterResp <- character(0)
    return(clusterResp)
  })
  names(haplotypes) <- sampleName
  tab <- table(unlist(haplotypes))
  overviewHap[names(tab),"nonchimera"] <- tab

  ######
  ## singleton $ representativ
  env <- environment()
  overviewHap$singelton <- T
  overviewHap$representatives <- F
  env$overviewHap$singelton <- overviewHap$singelton
  env$overviewHap$representatives <- overviewHap$representatives
  
  repfile <- clusterFilenames[,"RepresentativeFile"]
  haplotypes <- lapply(seq_along(repfile), function(i) {
  	sr1 <- readFasta(repfile[i])
  	vec <- do.call(rbind, strsplit(as.character(id(sr1)), "_"))
  	clusterResp <- vec[,1]
  	clusterSize <- as.integer(vec[,2])
  	env$overviewHap[clusterResp,"representatives"] <- T
  	clusterResp <- clusterResp[clusterSize>1]
  	env$overviewHap[clusterResp,"singelton"] <- F
  	return(NULL)
  })
  
  table(env$overviewHap$singelton)
  table(env$overviewHap$representatives)
  overviewHap$singelton <- env$overviewHap$singelton
  overviewHap$representatives <- env$overviewHap$representatives
  
  ######
  # INDEL - Homopolymere
	if(!is.null(referenceSequence)){
		idx <- as.character(id(sr1)) %in% as.character(overviewHap$HaplotypesName[!overviewHap$singelton])
		sr1 <- sr1[idx]
		aln1 <- pairwiseAlignment(sread(sr1), referenceSequence, type="global")
		idx <- (S4Vectors::lengths(deletion(aln1)) + S4Vectors::lengths(insertion(aln1))) == 0
		names(idx) <- as.character(id(sr1))
		overviewHap[names(idx),"indels"] <- !idx
		sr1 <- sr1[idx]
	}
  
  ###
  # only SNP variation - Final Haplotypes
  if(!is.null(snpList)){
	  snps <- lapply(as.integer(snpList[,"Pos"]), function(n) {
	    as.character(subseq(sread(sr1), start=n, width = 1))
	  })
	  snps <- apply(do.call(cbind, snps), 1, paste, collapse="")
	  names(snps) <- as.character(id(sr1))
	  
	  overviewHap[names(snps),"snps"] <- snps
  }else{
  	overviewHap[,"snps"] <- ""
  }
  return(overviewHap)
}
