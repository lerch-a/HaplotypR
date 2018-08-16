
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
  
  cat("\nbuilding contingency table...\n")
  freq <- integer(length(env$allHaplotypes))
  names(freq) <- names(env$allHaplotypes)
  contingencyTable <- lapply(contingencyList, function(l) {
    #idx <- match(l$UID, names(env$allHaplotypes))
    freq[as.character(l$UID)] <- l$Frequency
    return(freq)
  })
  print(nchar(env$allHaplotypes[1]))
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
                          minSampleCoverage=300, reportBackground=FALSE, defineBackground=NULL, verbose=FALSE, ...) {
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
  
  # remove low coverage sample
  cov <- colSums(x)
  pass_samplecov_filter_idx <- (cov >= minSampleCoverage)
  idx = pass_samplecov_filter_idx
  if (verbose) {
      print(paste(sum(pass_samplecov_filter_idx, na.rm=TRUE), "samples passed coverage filter."))
  } 
  if (all(!pass_samplecov_filter_idx)) {
    if (verbose) cat("\nNo samples passed sample coverage filter.")
  	x[,!idx] <- NA
  	#x <- x[1,, drop=F]
  	#rownames(x) <- NA
  	return(x)
  } else {
  	x[, !idx] <- NA 
  	if (verbose) {
  	    cat("\nafter sample coverage filter...\n")
  	    print(head(x))  
  	} 
  }

  # sample replicated
  minReplicate <- min(dim(x)[2], minReplicate)
  
  # remove haplotypes without reads
  x <- x[(rowSums(x > 0, na.rm=TRUE) > 0),, drop=F]
  if (verbose) {
      cat("\nafter removing haplotypes w/o reads...\n")
      print(head(x))  
  } 
  
  # selected and remove filtered haplotype
  idx <- rownames(x) %in% defineBackground
  background <- x[idx,, drop=F]
  x <- x[!idx,, drop=F]
  if (verbose) {
      cat("\nafter selecting pass-filtered haplotypes...\n")
      print(head(x))  
  } 
  
  # check for noise haplotype
  if (FALSE) { 
      minCov <- colSums(x) * detectability
      minCov[minCov < minHaplotypCoverage] <- minHaplotypCoverage
      if (dim(x)[1] > 0) {
      	noiseIdx <- (rowSums(t(t(x) / minCov) >= 1, na.rm=TRUE) < minReplicate) # only haplotypes present in min replicates
      	lowCnt <- colSums(x[noiseIdx,, drop=F])
      	x <- x[!noiseIdx,, drop=F]
      } else {
      	lowCnt <- 0
      }
      if (verbose) {
          cat("\nafter checking for noise haplotypes...\n")
          print(head(x))  
      } 
  }
  # add background to haplotyp counts
  if (reportBackground) {
    if (dim(x)[1] > 0)
        x <- rbind(Noise=lowCnt, background, x[order(rowSums(x)),, drop=F])
    else
        x <- rbind(Noise=lowCnt, background)
  }
  return(x)
}

# compute final table of called haplotypes 
createFinalHaplotypeTable <- function(outputDir, sampleTable, markerTab, snpLst, refSeq, postfix, 
                                      minHaplotypCoverage=2, minReplicate=3, detectability=0.01, minSampleCoverage=25, 
                                      include_seq=FALSE, verbose=FALSE, just_contingency_table=TRUE, flag_chimeras=FALSE,
                                      max_indel_thresholds=list(), return_full_haplotypes=FALSE) {
    #source("/Users/tfarrell/Tools/HaplotypR/R/clusterFunction.R")
    #source("/Users/tfarrell/Tools/HaplotypR/R/filterHaplotype.R")
  outFreqFiles <- file.path(outputDir, "haplotype_freq_files")
  dir.create(outFreqFiles)
  res <- lapply(rev(markerTab$MarkerID), function(marker) {
      if (verbose) cat(paste0("\nfor ", as.character(marker), ":\n"))
      max_indel_threshold = max_indel_thresholds[[as.character(marker)]]
    samTab <- sampleTable[sampleTable$MarkerID == marker,]
    potSNPLst <- snpLst[[marker]]
    prefix <- sub(".fastq.gz$", "", basename(as.character(samTab$ReadFile)))
    
    # Create frequency files and count table
    contingency_file = file.path(outputDir, sprintf("contingencyTable_%s%s.rds", marker, postfix[[marker]]))
    if (!file.exists(contingency_file)) { 
        if (verbose) cat("\ncomputing contingency table...\n")
        coverage_mat <- createContingencyTable(as.character(samTab$ReadFile), dereplicated=F, 
                                      inputFormat="fastq", outputDir=outFreqFiles,
                                      sampleNames=as.character(samTab$SampleID), replicatNames=NULL, 
                                      haplotypeUIDFile=NULL, include_seq=FALSE, verbose=verbose)
        if (verbose) cat("\nwriting contingency table...\n")
        saveRDS(coverage_mat, file=contingency_file)
    } else { 
        if (verbose) cat("\nreading in contingency table...\n")
        coverage_mat = readRDS(contingency_file) 
    }
    # compute uid coverage
    uid_cov = rowSums(coverage_mat[,sapply(colnames(coverage_mat), function (x) grepl("Solexa", x))])
    uid_cov = data.frame(uid_cov)
    colnames(uid_cov) = c("coverage")
    uid_cov["uid"] = rownames(uid_cov)
    
    # rename sequences fasta
    fnAllSeq <- file.path(outFreqFiles, sprintf("allSequences_%s%s.fasta", marker, postfix[[marker]]))
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
    
    if (!just_contingency_table) {
        overview_haplotype_table_file = file.path(outputDir, sprintf("HaplotypeOverviewTable_%s%s.txt", marker, postfix[[marker]]))
        if (!file.exists(overview_haplotype_table_file)) { 
            # create an overview table
            if (verbose) cat("\ncomputing haplotype overview table...\n")
            overviewHap <- createHaplotypOverviewTable(fnAllSeq, clusterFilenames, chimeraFilenames,
                                                       refSeq[[as.character(marker)]], potSNPLst, verbose=verbose, 
                                                       max_indel_threshold=max_indel_threshold)
            # label final haplotype as singleton, chimera, indel or {amplicon}-i for i in {1,...,N}, where N = # valid haplotypes
            overviewHap$FinalHaplotype <- factor(NA, levels = c("Singelton", "Chimera", "Indels", as.character(marker)))
            overviewHap[overviewHap$representatives, "FinalHaplotype"] <- as.character(marker)
            overviewHap[overviewHap$singelton, "FinalHaplotype"] <- "Singelton"
            overviewHap[!is.na(overviewHap$indels) & overviewHap$indels & overviewHap$representatives, "FinalHaplotype"] <- "Indels"
            overviewHap[!is.na(overviewHap$chimeraScore) & is.na(overviewHap$nonchimeraScore), "FinalHaplotype"] <- "Chimera"
            # cluster identical SNP patterns
            idx <- overviewHap$FinalHaplotype == as.character(marker)
            snps <- unique(overviewHap$snps[idx])
            if (length(snps) != 0) {
                names(snps) <- paste(marker, seq_along(snps), sep="-")
                levels(overviewHap$FinalHaplotype) <- c(levels(overviewHap$FinalHaplotype), names(snps))
                overviewHap$FinalHaplotype[idx] <- names(snps)[match(overviewHap$snps[idx], snps)]
            } 
            overviewHap <- as.data.frame(overviewHap)
            if (verbose) cat("\nwriting haplotype overview table...\n")
            write.table(overviewHap, file=overview_haplotype_table_file, sep="\t")
        } else { 
            if (verbose) cat("\nreading in haplotype overview table...\n")
            overviewHap = read.table(overview_haplotype_table_file)
        }
        if (verbose) {  
            cat("\noverview table:\n")
            print(head(overviewHap))
        } 
        suppressWarnings(suppressMessages(library(dplyr)))
        # compute filtered overview table, joined w/ uid coverage
        if (verbose) cat("\ncomputing filtered overview table...\n")
        overview_filtered = overviewHap[sapply(overviewHap["FinalHaplotype"], function (x) grepl(as.character(marker), x)), c("HaplotypesName","FinalHaplotype","snps")]
        colnames(overview_filtered) = c("uid","haplotypR_index","snps") 
        overview_filtered = merge(overview_filtered, uid_cov, by="uid")[,c("uid","haplotypR_index","coverage","snps")]
        # compute clustered overview table
        if (verbose) cat("\ncomputing clustered overview table...\n")
        overview_clustered = group_by(overview_filtered, haplotypR_index)
        overview_clustered = data.frame(summarize(overview_clustered, total_coverage=sum(coverage), snps=snps[which.max(coverage)]))
        if (return_full_haplotypes) { 
            # convert snps to full haplotype
            if (verbose) cat("\nconverting SNPs to full haplotypes...\n")
            ref = as.character(refSeq[[as.character(marker)]])
            snp_poss = snpLst[[as.character(marker)]]$Pos
            if (verbose) { 
                print(snp_poss)
                print(ref)
            } 
            
            for (r in rownames(overview_clustered)) { 
                h = ref
                s = as.character(overview_clustered[r, 'snps'])
                if (!is.na(s)) { 
                    for (i in seq_along(snp_poss)) { substr(h, snp_poss[i], snp_poss[i]) = substr(s, i, i) }
                } 
                overview_clustered[r, 'haplotype'] = h
            } 
        } 
        # write to file 
        write.table(overview_clustered, file=paste0("haplotypR.", as.character(marker), ".haplotype.index.tsv"), 
                    sep='\t', quote=FALSE)
        detach("package:dplyr", unload=TRUE)
        rm(list=c("overview_filtered", "uid_cov", "coverage_mat"))
        
        raw_final_haplotype_table_file = file.path(outputDir, sprintf("rawHaplotypeTable_%s%s.rds", marker, postfix[[marker]]))
        if (!file.exists(raw_final_haplotype_table_file)) { 
            if (verbose) cat("\nreading in representative files for final haplotype table...\n")
            hab <- rep(0, dim(overviewHap)[1])
            names(hab) <- rownames(overviewHap)
            repfile <- clusterFilenames[,"RepresentativeFile"]
            haplotypesSample <- lapply(seq_along(repfile), function(i) {
                if (verbose) cat(paste0("retrieving ", repfile[i], "...\n"))
              sr1 <- readFasta(repfile[i])
              vec <- do.call(rbind, strsplit(as.character(id(sr1)), "_"))
              clusterResp <- vec[,1]
              clusterSize <- as.integer(vec[,2])
              hab[clusterResp] <- clusterSize
              hab
            })
            if (verbose) cat("\nbuilding raw final haplotype table...\n")
            haplotypesSample <- do.call(cbind, haplotypesSample)
            haplotypesSample <- haplotypesSample[rowSums(haplotypesSample) > 1,]
            colnames(haplotypesSample) <- sub(".representatives.fasta", "", basename(repfile))
            overviewHap <- overviewHap[rownames(haplotypesSample),]
            rownames(haplotypesSample) <- overviewHap[rownames(haplotypesSample), "FinalHaplotype"]
            if (verbose) print(head(haplotypesSample))
            if (verbose) cat("\nwriting raw final haplotype table...\n")
            saveRDS(cbind(HaplotypNames=rownames(haplotypesSample), haplotypesSample), 
                    file=raw_final_haplotype_table_file)
        }
        reclustered_final_haplotype_table_file = file.path(outputDir, sprintf("reclusteredHaplotypeTable_%s%s.rds", marker, postfix[[marker]]))
        if(!file.exists(reclustered_final_haplotype_table_file)) {
            if (verbose) cat("\nretrieving raw final haplotype table...\n")
            haplotypesSample = readRDS(raw_final_haplotype_table_file)
            haplotype_names = rownames(haplotypesSample)
            haplotypesSample = apply(haplotypesSample, 2, as.integer)
            rownames(haplotypesSample) = haplotype_names
            if (verbose) cat("\nreclustering final haplotype table...\n")
            haplotypesSample <- reclusterHaplotypesTable(haplotypesSample)
            saveRDS(cbind(HaplotypNames=rownames(haplotypesSample),haplotypesSample), 
                    file=reclustered_final_haplotype_table_file)
        } 
        called_final_haplotype_table_file = file.path(outputDir, sprintf("finalCalledHaplotypeTable_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.rds", 
                                                                         minHaplotypCoverage, minSampleCoverage, minReplicate, detectability, 
                                                                         as.character(marker), postfix[[marker]]))
        if (!file.exists(called_final_haplotype_table_file)) {
            if (verbose) cat("\nretrieving reclustered final haplotype table...\n")
            haplotypesSample = readRDS(reclustered_final_haplotype_table_file) 
            haplotype_names = rownames(haplotypesSample)
            haplotypesSample = apply(haplotypesSample, 2, as.integer)
            rownames(haplotypesSample) = haplotype_names
            #if (verbose) print(head(haplotypesSample))
            if (verbose) cat("\ncalling final haplotype table...\n")
            haplotypesSample <- callHaplotype(haplotypesSample, minHaplotypCoverage=minHaplotypCoverage, minReplicate=minReplicate, 
                                              detectability=detectability, minSampleCoverage=1)
            saveRDS(cbind(HaplotypNames=rownames(haplotypesSample), haplotypesSample), 
                    file=called_final_haplotype_table_file)
        } else { 
            if (verbose) cat("\nretrieving called final haplotype table...\n")
            haplotypesSample = readRDS(called_final_haplotype_table_file) 
            haplotype_names = rownames(haplotypesSample)
            haplotypesSample = apply(haplotypesSample, 2, as.integer)
            rownames(haplotypesSample) = haplotype_names
            if (verbose) { 
                cat("\ncalled final haplotypes table:\n")
                print(head(haplotypesSample))
            }
        }
        if (verbose) cat("\nchecking final haplotypes for # replicates and chimeras...\n")
        idx <- split(1:dim(haplotypesSample)[2], samTab$SampleName)
        idx = idx[[names(idx)[1]]]
        if (verbose) cat(paste0("\nrecalling haplotypes...\n"))
        tab <- callHaplotype(haplotypesSample[,idx, drop=F], minHaplotypCoverage=minHaplotypCoverage, 
                             minReplicate=length(idx), detectability=detectability, minSampleCoverage=minSampleCoverage, 
                             reportBackground=FALSE, verbose=TRUE)
        if (verbose) { 
            cat("recalled haplotypes...\n")
            print(head(tab))
        }
        # if no coverage for any samples, output empty finalHaplotypeList
        if (all(colSums(tab) == 0)) { 
            result = setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Solexa_ID", "haplotype_index", "coverage"))
            write.table(result, file=file.path(outputDir, sprintf("finalHaplotypeList_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt", 
                                                                  minHaplotypCoverage, minSampleCoverage, minReplicate, detectability, 
                                                                  as.character(marker), postfix[[marker]])), 
                        sep="\t", row.names=F, col.names=T, quote=FALSE)
            return(result)
        }
        if (verbose) cat(paste0("\nmelting haplotypes table...\n"))
        tab = melt.matrix(tab)
        colnames(tab) = c("haplotypR_index","sample_id","coverage")
        tab = tab[!is.na(tab["coverage"]),]
        tab = tab[tab["coverage"] > 0,]
        
        if (flag_chimeras) { 
          if (verbose) cat("\nchecking for chimeras...\n")
          if (length(idx)) {
            tmpTab <- as.data.frame(tapply(tab$Reads, tab$Haplotype, sum))
            tmpTab$Haplotype <- rownames(tmpTab)
            colnames(tmpTab) <- c("Reads","Haplotype")
            hIdx <- grep(marker, tmpTab$Haplotype)
          } else {
            hIdx <- grep(marker, tab$Haplotype)
          }
          chim <- NULL
          if(length(hIdx) > 2){
            chim <- flagChimera(tmpTab[hIdx,], overviewHap)
          }
          rm(tmpTab)
          tab$FlagChimera <- tab$Haplotype %in% chim
        }
        if (verbose) {
            cat("\nfinal haplotypes table:\n")
            print(head(tab))
        }
        
        # add Solexa_ID to final list 
        tab$Solexa_ID = sapply(tab$sample_id, function(x) substr(x, 0, 13))
        # join full haplotypes to final list
        if ("haplotype" %in% colnames(overview_clustered)) { 
            full_tab = merge(tab[,c("Solexa_ID","coverage","haplotypR_index")], 
                             overview_clustered[,c("haplotypR_index","haplotype")], 
                             by="haplotypR_index", all.x=TRUE)
        } else { 
            full_tab = merge(tab[,c("Solexa_ID","coverage","haplotypR_index")], 
                             overview_clustered[,c("haplotypR_index","snps")], 
                             by="haplotypR_index", all.x=TRUE)    
        }
        # write to file
        cat("\nwriting final haplotype table...\n")
        write.table(full_tab, file=file.path(outputDir, sprintf("finalHaplotypeList_Hcov%.0f_Scov%.0f_occ%i_sens%.4f_%s%s.txt", 
                                                           minHaplotypCoverage, minSampleCoverage, minReplicate, detectability, 
                                                           as.character(marker), postfix[[marker]])), 
                    sep="\t", row.names=F, col.names=T, quote=FALSE)
        rownames(tab) <- NULL
        return(tab)
        } else { 
            return(contingency_file)    
        }
    })
  names(res) <- markerTab$MarkerID
  return(res)
}