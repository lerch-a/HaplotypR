
createProject <- function(){
  
  list(outDir=outputDir,
       sampleFile=sampleFile,
       primerFile=primerFile,
       marker=primer$MarkerID,
       refSeq=primer$ReferenceSequence,
       projects=prj
  )
}

# createProject <- function(outputDir, primerFile, sampleFile=NULL){
#   
#   primer <- read.delim(primerFile, stringsAsFactors=F)
#   fn <- list.files(outputDir, "processedReadSummary")
#   if(length(fn)>0){
#     prj <- do.call(rbind, strsplit(fn, "_bind"))[,2]
#     prj <- sub(".txt", "", prj)
#     names(prj) <-  paste("_bind", prj, sep="")
#     prj <- strsplit(prj, "_")
#     if(length(idx <- grep("full", prj))>0)
#       prj[[idx]] <- c(NA, NA)
#     prj <- as.data.frame(do.call(rbind, prj))
#     colnames(prj) <- c("numNtFwd","numNtRev")
#     prj$filename <- fn
#   }else{
#     prj <- NULL
#   }
# 
#   list(outDir=outputDir,
#        sampleFile=sampleFile,
#        primerFile=primerFile,
#        marker=primer$MarkerID,
#        refSeq=primer$ReferenceSequence,
#        projects=prj
#   )
# }
