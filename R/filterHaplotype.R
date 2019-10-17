

checkChimeras <- function(representativesFile, method="vsearch", progressReport=message){
  
  # sort representatives 
  sortfile <- unlist(lapply(seq_along(representativesFile), function(i){
    #vsearch --sortbysize R-1_Rep1.representatives.fasta --output R-1_Rep1.representatives_sortbysize.fasta
    sr1 <- readFasta(representativesFile[i])
    clusterSize <- as.integer(do.call(rbind, strsplit(as.character(id(sr1)), "_"))[,2])
    sr1 <- sr1[order(clusterSize, decreasing=T)]
    sortfile <- sub(".fasta", "_sort.fasta", representativesFile[i])
    writeFasta(ShortRead(sread(sr1), BStringSet(sub("_", ";size=", as.character(id(sr1))))), sortfile)
    return(sortfile)
  }))
  
  nonchimerafile <- sub("_sort.fasta", "_nonchimera.fasta", sortfile)
  chimerafile <- sub("_sort.fasta", "_chimera.fasta", sortfile)
  borderfile <- sub("_sort.fasta", "_bordchimera.fasta", sortfile)
  resfile <- sub("_sort.fasta", "_chimeraResults.txt", sortfile)
  # vsearch --uchime_denovo R-1_Rep1.representatives_nonchimeras.fasta --nonchimeras R-1_Rep1.representatives_nonchimeras2.fasta
  syscall <- paste("--uchime_denovo ", sortfile,"--mindiffs", 3, "--minh", 0.2, "--nonchimeras", nonchimerafile, "--chimeras", chimerafile, "--borderline", borderfile, "--uchimeout", resfile, sep=" ")
  
  # check and set progress report function
  if(!is.function(progressReport))
    progressReport <- message
  
  require(Rvsearch)
  #Rvsearch:::.vsearchBin(args=syscall)
  lapply(seq_along(syscall), function(i){ 
    msg <- paste("Processing file", basename(representativesFile[i]), "...", sep=" ")
    progressReport(detail=msg, value=i)
    Rvsearch:::.vsearchBin(args=syscall[i])
  })
  
  return(cbind(NonchimeraFile=nonchimerafile, ChimeraFile=chimerafile, BorderchimeraFile=borderfile, ChimeraResultsFile=resfile))
}


createHaplotypOverviewTable <- function(allHaplotypesFilenames, clusterFilenames, chimeraFilenames, referenceSequence,
                                        snpSet, maxDel = 9L, maxIns = 9L){
  
  sr1 <- readFasta(allHaplotypesFilenames)
  overviewHap <- data.frame(HaplotypesName=as.character(id(sr1)))
  rownames(overviewHap) <- overviewHap$HaplotypesName

  ######
  ## chimera type
  resfile <- chimeraFilenames[,"ChimeraResultsFile"]
  #sampleName <- sub(".representatives_chimeraResults.fasta", "", basename(resfile))
  haplotypes <- lapply(seq_along(resfile), function(i){
  	res <- read.delim(resfile[i], header = F, col.names=paste("V", rep(1:18), sep="")) # colnames that empty imput does not give error
  	chimScore <- res[,1]
  	vec <- do.call(rbind, strsplit(as.character(res[,2]), ";size="))
  	clusterResp <- vec[,1]
  	clusterSize <- as.integer(vec[,2])
  	clusterResp <- clusterResp[clusterSize>1]
  	chimScore <- chimScore[clusterSize>1]
  	if (length(clusterResp)==0){
  		return(list(Chimera=character(0), NonChimera=character(0)))
  	}
  	return(list(Chimera=clusterResp[chimScore>0], NonChimera=clusterResp[chimScore==0]))
  })
  #names(haplotypes) <- sampleName
  chim <- table(unlist(lapply(haplotypes, "[[", 1)))
  nochim <- table(unlist(lapply(haplotypes, "[[", 2)))
  
  overviewHap[names(chim),"chimeraScore"] <- chim
  overviewHap[names(nochim),"nonchimeraScore"] <- nochim
  
  chimerafile <- chimeraFilenames[,"ChimeraFile"]
  sampleName <- sub(".representatives_chimera.fasta", "", basename(chimerafile))
  #sampleName <- do.call(rbind, strsplit(basename(chimerafile), "\\."))[,1]
  haplotypes <- lapply(seq_along(chimerafile), function(i){
    sr1 <- readFasta(chimerafile[i])
    vec <- do.call(rbind, strsplit(as.character(id(sr1)), ";size="))
    clusterResp <- vec[,1]
    clusterSize <- as.integer(vec[,2])
    clusterResp <- clusterResp[clusterSize>1]
    if(length(clusterResp)<0)
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
    if(length(clusterResp)<0)
      clusterResp <- character(0)
    return(clusterResp)
  })
  names(haplotypes) <- sampleName
  tab <- table(unlist(haplotypes))
  overviewHap[names(tab),"bordchimera"] <- tab
  
  
  nonchimerafile <- chimeraFilenames[,"NonchimeraFile"]
  haplotypes <- lapply(seq_along(nonchimerafile), function(i){
    sr1 <- readFasta(nonchimerafile[i])
    vec <- do.call(rbind, strsplit(as.character(id(sr1)), ";size="))
    clusterResp <- vec[,1]
    clusterSize <- as.integer(vec[,2])
    clusterResp <- clusterResp[clusterSize>1]
    if(length(clusterResp)<0)
      clusterResp <- character(0)
    return(clusterResp)
  })
  names(haplotypes) <- sampleName
  tab <- table(unlist(haplotypes))
  overviewHap[names(tab),"nonchimera"] <- tab

  
  ######
  ## singleton $ representativ
  overviewHap$singelton <- T
  overviewHap$representatives <- F
  
  repfile <- clusterFilenames[,"RepresentativeFile"]
  for (fn in repfile){
    sr1 <- readFasta(fn)
    vec <- do.call(rbind, strsplit(as.character(id(sr1)), "_"))
    at <- match(vec[,1], overviewHap$HaplotypesName)
    clusterSize <- as.integer(vec[,2])
    overviewHap$representatives[at] <- T
    overviewHap$singelton[at[clusterSize>1]] <- F
  }
  
  #  Indels
  overviewHap$indels <- F
  hap_set <- as.character(with(overviewHap, HaplotypesName[representatives & ! singelton]))
	haps <- readDNAStringSet(allHaplotypesFilenames)[hap_set]
	aln1 <- pairwiseAlignment(haps, referenceSequence, type="global")
	del_width <- as.list(deletion(aln1)) %>% purrr::map_int(~ sum(width(.)))
	ins_width <- as.list(insertion(aln1)) %>% purrr::map_int(~ sum(width(.)))
	indels_at <- which(del_width > maxDel | ins_width > maxIns)
	overviewHap$indels[match(hap_set, overviewHap$HaplotypesName)[indels_at]] <- T
	
	# SNPs
	overviewHap$snps <- rep(NA_character_, nrow(overviewHap))
  if (length(indels_at) > 0){
    hap_set <- hap_set[-indels_at]
    haps <- haps[hap_set]
  }
	aln2 <- pairwiseAlignment(haps, referenceSequence, type="global")
  snp_df <-
    purrr::map_df(seq_along(hap_set), function(i) {
      mismatchSummary(aln2[i])$subject %>% 
        {`if`('Pattern' %in% colnames(.), ., dplyr::mutate(., Pattern = character()))} %>% 
        dplyr::select(Pos=SubjectPosition, Ref=Subject, Alt=Pattern) %>% 
        dplyr::mutate_if(is.factor, as.character) %>% 
        dplyr::right_join(snpSet, c('Pos', 'Ref')) %>% 
        dplyr::summarise(snp = dplyr::if_else(is.na(Alt), Ref, Alt) %>% paste0(collapse = '')) %>%
        dplyr::mutate(uid = hap_set[i])
    })
  
  overviewHap$snps[match(snp_df$uid, overviewHap$HaplotypesName)] <- snp_df$snp
  return(overviewHap)
}
