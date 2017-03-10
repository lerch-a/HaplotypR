
plotMismatchFrequencies <- function(control, sample, filename, ...){
  
  seqErrFwdC <- do.call(cbind, lapply(control, function(l){ 
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  
  seqErrFwd <- do.call(cbind, lapply(sample, function(l){ 
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  
  png(filename, width=1500 , height=600)
  matplot(seqErrFwd, type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1), 
          ylab="Mismatch Rate", xlab="nt position", main="CSP Fragment", cex.axis=1.5, cex.lab=1.5)
  matplot(seqErrFwdC, type="p", pch=16, cex=0.4, col="#FF000088", add=T)
  abline(v=190.5, lty=1, col="black")
  potSNP <- rowSums(seqErrFwd>0.5)>=2
  potSNP <- seq_along(potSNP)[potSNP]
  abline(v=potSNP, lty=2, col="grey")
  abline(h=0.5, lty=3, col="red")
  dev.off()
  
}


plotHaplotypes <- function(allelCounts, sampleLable=colnames(allelCounts), title, minCoverage=3, maxSensitivity=1/1000){
  
  allelCounts <- callHaplotype(allelCounts, sensitivity=maxSensitivity, minCoverage=minCoverage)
  
  totalCoverage <- colSums(allelCounts)
   
  # # remove haplotypes without reads
  # idx <- rowSums(allelCounts>0) > 0
  # allelCounts <- allelCounts[idx,,drop=F]
  
  #
  # idx <- rowSums(allelCounts>=minCounts) > 0
  # lowCnt <- colSums(allelCounts[!idx,,drop=F])
  # allelCounts <- allelCounts[idx,,drop=F]
  # allelCounts <- allelCounts[order(rowSums(allelCounts), decreasing=F),,drop=F]
  
  # if(dim(allelCounts)[1]>0)
  #   allelCounts <- rbind(Noise=lowCnt, allelCounts)
  # else
  #   allelCounts <- t(data.frame(Noise=lowCnt))
  
  # nIdx <- grep("Noise", rownames(allelCounts))
  # rownames(allelCounts)[nIdx] <- sprintf("%s <%ix", rownames(allelCounts)[nIdx], minCounts)
  
  # colnames(allelCounts) <- sampleLable
  # #sort sample by lable
  # allelCounts <- allelCounts[,order(colnames(allelCounts)),drop=F]
  
  xlab <- colnames(allelCounts)
  ylab <- c(rownames(allelCounts), "Coverage")
  
  plot(NA, ylim=c(0.5, dim(allelCounts)[1]+1.5), xlim=c(0.5,dim(allelCounts)[2]+0.5),
       type="n", bty="n", axes = FALSE, 
       xlab="", ylab="", main=title)
  
  axis(side=2, at=1:(dim(allelCounts)[1]+1), labels=ylab, tick=F, las=2, cex.axis=0.8)
  axis(side=1, at=1:dim(allelCounts)[2], labels=xlab, tick=F, hadj=0.7, padj=0.5, las=2, cex.axis=0.8)
  
  if(dim(allelCounts)[1]>0){ # plot only if haplotype exist
    invisible(lapply(1:dim(allelCounts)[2], function(i){
      s <- sum(allelCounts[,i])
      if(s>0){ # check if any reads in cluster
        perc <- allelCounts[,i]
        perc[perc==0] <- NA
        symbols(rep(i,dim(allelCounts)[1]), 1:dim(allelCounts)[1], 
                circles=(perc/s)*(dim(allelCounts)[2]/20), 
                inches=F, add=T, fg="#FF0000FF", bg="#FF0000FF")
        
        text(i, 1:dim(allelCounts)[1], labels=sprintf("%.1f%%", perc/s*100), cex=0.7, offset=0)
      }
      text(i, dim(allelCounts)[1]+1, labels=sprintf("%sx", format(s, big.mark="'")), cex=0.7, offset=0)
      return(NULL)
    }))
  }
}
