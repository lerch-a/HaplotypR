

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
                                        snpSet, maxDel = 0L, maxIns = 0L){
  
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
	del_width <- vapply(as.list(deletion(aln1)), function(x) sum(width(x)), numeric(1))
	ins_width <- vapply(as.list(insertion(aln1)), function(x) sum(width(x)), numeric(1))
	indels_at <- which(del_width > maxDel | ins_width > maxIns)
	overviewHap$indels[match(hap_set, overviewHap$HaplotypesName)[indels_at]] <- T
	
	# SNPs
	overviewHap$snps <- NA_character_
  if (length(indels_at) > 0){
    hap_set <- hap_set[-indels_at]
    haps <- haps[hap_set]
  }
	aln2 <- pairwiseAlignment(haps, referenceSequence, type="global")
  snp_df <-
    lapply(seq_along(hap_set), function(i) {
      df1 <- mismatchSummary(aln2[i])$subject
      df1 <- `if`('Pattern' %in% colnames(df1), df1, { df1$Pattern <- character(); df1 })
      df1$Subject <- as.character(df1$Subject)
      df1$Pattern <- as.character(df1$Pattern)
      df1 <- data.frame(Pos = df1$SubjectPosition, Ref = df1$Subject, Alt = df1$Pattern)
      df1 <- merge.data.frame(df1, snpSet, by = c('Pos', 'Ref'), all.y = TRUE, sort = FALSE)
      df1$snp <- ifelse(is.na(df$Alt), df$Ref, df$Alt)
      df2 <- data.frame(uid = hap_set[i], snp = paste0(df1$snp, collapse = ''))
      return(df2)
    })
  snp_df <- do.call(rbind.data.frame, snp_df)
  
  overviewHap$snps[match(snp_df$uid, overviewHap$HaplotypesName)] <- snp_df$snp
  return(overviewHap)
}
