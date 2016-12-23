

reclusterHaplotypesTable <- function(x, f=NULL){
  hapName <- if(is.null(f)){
    unique(rownames(haplotypesSample))
  }else{
    stopifnot(length(f) == dim(x)[1])
    f
  }
    
  x <- do.call(rbind, lapply(hapName, function(name){
    colSums(x[rownames(x) %in% name,, drop=F])
  }))
  rownames(x) <- hapName
  return(x)
}

