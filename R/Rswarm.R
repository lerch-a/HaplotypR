
## Example
# syscall <- paste("-f -b 3 -w", repfile, "-l", logfile, "-s", statsfile, "-i", structfile, "-o", outfile, frqfile, sep=" ")
# swarm(syscall)

## The main wrapper around swarm
swarm <- function(fastaFile, outfilePrefix=sub("\\.fasta", "", fastaFile), 
                  boundary=3, differences=1, fastidious=T, noOTUbreaking=F,
                  force=F, internalStructure=T, log=T, outputFile=T, statisticsFile=T, seeds=T, 
                  mothurFormat=F, usearchAbundanceStyle=F, appendAbundance=NULL,
                  matchReward=5, mismatchPenalty=4, gapOpeningPenalty=12, gapExtensionPenalty=4,
                  threads=1, ceiling=NULL, bloomBits=16)
{
  args <- ""
  
  ##    General options:
  #     -t, --threads INTEGER 1-256               number of threads to use (1)
  if(threads %in% 1:256)
    args <- paste(args, sprintf("-t %d", threads))
  
  
  ##   Clustering options:
  #     -b, --boundary INTEGER              min mass of large OTU for fastidious (3)
  if(is.integer(boundary))
    args <- paste(args, sprintf("-b %d", boundary))
  
  #     -c, --ceiling INTEGER               max memory in MB used for fastidious
  if(!is.null(ceiling) & is.integer(ceiling))
    args <- paste(args, sprintf("-c %d", ceiling))
  
  #     -d, --differences INTEGER 0-256           resolution (1)
  if(differences %in% 0:256)
    args <- paste(args, sprintf("-d %d", differences))
  
  #     -f, --fastidious                    link nearby low-abundance swarms
  if(fastidious)
    args <- paste(args, "-f")
  
  #     -n, --no-otu-breaking               never break OTUs
  if(noOTUbreaking)
    args <- paste(args, "-n")
  
  #     -y, --bloom-bits INTEGER 2-64            bits used per Bloom filter entry (16)
  if(bloomBits %in% 2:64)
    args <- paste(args, sprintf("-y %d", bloomBits))
  
  
  ##    Input/output options:
  #     -a, --append-abundance INTEGER      value to use when abundance is missing
  if(!is.null(appendAbundance) & is.integer(appendAbundance))
    args <- paste(args, sprintf("-a %d", appendAbundance))
    
  #     -i, --internal-structure FILENAME   write internal swarm structure to file
  if(internalStructure){
    structfile <- paste(outfilePrefix, ".struct.txt", sep="")
    args <- paste(args, sprintf("-i %s", structfile))
  }
  
  #     -l, --log FILENAME                  log to file, not to stderr
  if(log){
    logfile <- paste(outfilePrefix, ".log.txt", sep="")
    args <- paste(args, sprintf("-l %s", logfile))
  }
    
  #     -o, --output-file FILENAME          output result filename (stdout)
  if(outputFile){
    outfile <- paste(outfilePrefix, ".swarms.txt", sep="")
    args <- paste(args, sprintf("-o %s", outfile))
  }
  
  #     -r, --mothur                        output in mothur list file format
  if(mothurFormat)
    args <- paste(args, "-r")
  
  #     -s, --statistics-file FILENAME      dump OTU statistics to file
  if(statisticsFile){
    statsfile <- paste(outfilePrefix, ".stats.txt", sep="")
    args <- paste(args, sprintf("-s %s", statsfile))
  }
  
  #     -u, --uclust-file FILENAME          output in UCLUST-like format to file
  # TODO
  
  #     -w, --seeds FILENAME                write seed seqs with abundances to FASTA
  if(seeds){
    repfile <- paste(outfilePrefix, ".representatives.fasta", sep="")
    args <- paste(args, sprintf("-w %s", repfile))
  }
  
  #     -z, --usearch-abundance             abundance annotation in usearch style
  if(usearchAbundanceStyle)
    args <- paste(args, "-z")
  
  
  ##    Pairwise alignment advanced options:
  #     -m, --match-reward INTEGER          reward for nucleotide match (5)
  if(is.integer(matchReward))
    args <- paste(args, sprintf("-m %d", matchReward))
  
  #     -p, --mismatch-penalty INTEGER      penalty for nucleotide mismatch (4)
  if(is.integer(mismatchPenalty))
    args <- paste(args, sprintf("-p %d", mismatchPenalty))
                  
  #     -g, --gap-opening-penalty INTEGER   gap open penalty (12)
  if(is.integer(gapOpeningPenalty))
    args <- paste(args, sprintf("-g %d", gapOpeningPenalty))
                  
  #     -e, --gap-extension-penalty INTEGER gap extension penalty (4)  
  if(is.integer(gapExtensionPenalty))
    args <- paste(args, sprintf("-e %d", gapExtensionPenalty))

  
  if(!is.character(fastaFile) || file.exists(fastaFile)){
    stop("Argument 'fastaFile' has to be a character vector of filenames.")
  }
  args <- paste(args, fastaFile)

  return(invisible(.swarmBin(args)))
}


## Helper function that return a description of the intended usage for swarm
swarm_usage <- function()
  print(.swarmBin(args="--help"))


## Helper function that return the version of swarm
swarm_version <- function(){
  print(.swarmBin(args="--version"))
}


## A helper function to call the swarm binaries with additional arguments.
.swarmBin <- function(args="")
{
  if(is.null(args) || args=="")
    stop("The swarm binaries need to be called with additional arguments")
  args <- gsub("^ *| *$", "", args)
  
  swarmExec <- list.files(system.file(file.path('swarm', 'bin'), package = 'HaplotypR'),
                          pattern = 'swarm', full.names = TRUE)[1]
  call <- paste(shQuote(swarmExec), args)
  
  return(system(call, intern=TRUE))
}

## The direct binary call function
.swarmExecute <- function(callstr, ...){
  
  swarmExec <- list.files(system.file(file.path('swarm', 'bin'), package = 'HaplotypR'),
                          pattern = 'swarm', full.names = TRUE)[1]
  call <- file.path(shQuote(swarmExec), callstr)
  
  return(system(call, ...))
}

