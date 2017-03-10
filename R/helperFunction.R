

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




