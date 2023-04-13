
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

reclusterHaplotypesTable <- function(haplotypeTable, haplotypCluster=NULL){
  hapName <- if(is.null(haplotypCluster)){
    unique(rownames(haplotypeTable))
  }else{
    stopifnot(length(haplotypCluster) == dim(haplotypeTable)[1])
    rownames(haplotypeTable) <- haplotypCluster
    unique(haplotypCluster)
  }
    
  haplotypeTable <- do.call(rbind, lapply(hapName, function(name){
    colSums(haplotypeTable[rownames(haplotypeTable) %in% name,, drop=F])
  }))
  rownames(haplotypeTable) <- hapName
  return(haplotypeTable)
}

compareDNAString <- function(query, subject){
	query <- as.character(query)
	subject <- as.character(subject)
	bases <- c("A","C","G","T")
	mat <- consensusMatrix(query)[bases,]
	refSeq <- strsplit(subject, "")[[1]]
	refIdx <- match(refSeq, bases)
	
	mismatchSummary <- lapply(seq_along(refSeq), function(i){
		cov <- sum(mat[,i])
		#mat[refIdx[i],i]
		mm <- mat[-refIdx[i],i]
		data.frame(SubjectPosition=i, Subject=refSeq[i], Pattern=names(mm), Count=mm, Probability=mm/cov)
	})
	mismatchSummary <- do.call(rbind, mismatchSummary)
	rownames(mismatchSummary) <- NULL
	mismatchSummary <- mismatchSummary[mismatchSummary$Count>0,]
	return(mismatchSummary)
}

findChimeras <- function(snps){
	
	alnLst <- unlist(lapply(1:(length(snps)-1), function(i){
		unlist(lapply((i+1):length(snps), function(j){
			aln <- compareStrings(snps[j], snps[i])
			names(aln) <- paste(i, j)
			return(aln)
		}) )
	}))
	
	res <- do.call(rbind, lapply(3:length(snps), function(j){
		sel <- alnLst[grep(sprintf(". %i", j), names(alnLst))]
		pos = gregexpr('\\?', sel)
		names(pos) <- names(sel)
		len <- lengths(pos)
		comb <- combn(names(sel), 2)
		do.call(rbind, lapply(1:dim(comb)[2], function(i){
			if((pos[[comb[,i][1]]][1] > pos[[comb[,i][2]]][1]) & (pos[[comb[,i][1]]][1] >= pos[[comb[,i][2]]][len[comb[,i][2]]])){
				src <- names(snps)[as.integer(substring(comb[,i],1,1))]
				chim <- names(snps)[as.integer(substring(comb[1,i],3,3))]
				return(c(chim, src))
			}else if((pos[[comb[,i][1]]][1] < pos[[comb[,i][2]]][1]) & (pos[[comb[,i][1]]][len[comb[,i][1]]] <= pos[[comb[,i][2]]][1])){
				src <- names(snps)[as.integer(substring(comb[,i],1,1))]
				chim <- names(snps)[as.integer(substring(comb[1,i],3,3))]
				return(c(chim, src))
			}else return(NULL)
		}))
	}))
	if(length(res)>0)
		colnames(res) <- c("chimera", "src1", "src2")
	return(res)
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


getReferenceSequence <- function(genome, forwardPrimer, reversePrimer){
  refSeq <- apply(cbind(Forward=forwardPrimer, Reverse=reversePrimer), 1, function(rr){
    Fprobe <- rr["Forward"]
    Rprobe <- rr["Reverse"]
    # find amplicon
    lstViews <- lapply(genome, function(sub){
      matchProbePair(as.character(Fprobe), as.character(Rprobe), sub, verbose=F) # 1 "theoretical amplicon"
    })
    lstViews <- lstViews[lengths(lstViews)>0]
    if(length(lstViews)!=1) # check if single match
      return("")
    seq <- unlist(lapply(lstViews, function(vw){ 
      toString(vw)
    }))
    # check strand
    if(substr(seq, 1, nchar(as.character(Fprobe))) != as.character(Fprobe))
      seq <- toString(reverseComplement(DNAString(seq)))
    # remove primer sequence
    seq <- subseq(seq, nchar(as.character(Fprobe))+1, width(seq))
    seq <- subseq(seq, 1, width(seq)-nchar(as.character(Rprobe)))
    return(seq)
  })
}

getReferenceLocation <- function(genome, forwardPrimer, reversePrimer, includPrimer=F){
  names(genome) <- do.call(rbind, strsplit(names(genome)," "))[,1]
  refLocation <- apply(cbind(Forward=forwardPrimer, Reverse=reversePrimer), 1, function(rr){
    Fprobe <- rr["Forward"]
    Rprobe <- rr["Reverse"]
    # find amplicon
    lstViews <- lapply(genome, function(sub){
      matchProbePair(as.character(Fprobe), as.character(Rprobe), sub, verbose=F) # 1 "theoretical amplicon"
    })
    lstViews <- lstViews[lengths(lstViews)>0]
    if(length(lstViews)!=1) # check if single match
      return("")
    seq <- toString(lstViews[[1]])
    # check strand
    isFwdStrand <- substr(seq, 1, nchar(as.character(Fprobe))) == as.character(Fprobe)
    gr <- GRanges(seqnames=names(lstViews)[1], ranges=ranges(lstViews[[1]]), strand = ifelse(isFwdStrand,"+","-"))
    # # remove primer sequence
    if(!includPrimer){
      if(isFwdStrand){
        # seq <- subseq(seq, nchar(as.character(Fprobe))+1, width(seq))
        start(gr) <- start(gr) + nchar(as.character(Fprobe))
        end(gr) <- end(gr) - nchar(as.character(Rprobe))
        subseq(genome[seqnames(gr)], start(gr), end(gr))
      }else{
        start(gr) <- start(gr) + nchar(as.character(Rprobe))
        end(gr) <- end(gr) - nchar(as.character(Fprobe))
        subseq(genome[seqnames(gr)], start(gr), end(gr))
      }
    }
    return(gr)
  })
  return(unlist(as(refLocation,"GRangesList")))
}

assignStrain <- function(haplotypes, refFn){
  refSeq <- readDNAStringSet(refFn)
  #refSeq <- narrow(refSeq, 1, width(haplotypes)[1])
  refName <- unlist(lapply(haplotypes, function(hap){
    # hap <- haplotypes[[1]]
    paste(names(refSeq)[hap==refSeq], collapse=";")
  }))
  return(refName)
}

assignSNP <- function(haplotypes, snpTab){
  snpTab <- snpTab[snpTab$ALT!="",]
  if(nrow(snpTab)>0){
    mutLst <- do.call(rbind, lapply(1:nrow(snpTab), function(ii){
      message(ii)
      pos <- snpTab$nt_position_coding_strand[ii]
      snp <- subseq(haplotypes, pos, width=1)
      ifelse(snpTab$ALT[ii]==as.character(snp), snpTab$Mutation[ii], paste0("p",pos))
    }))
    mutLst <- apply(mutLst, 2, function(hh){
      hh <- hh[!is.na(hh)]
      paste0(hh, collapse=",")
    })
    return(mutLst)
  }else return(NULL)
}

