#
# adapted from https://github.com/lerch-a/HaplotypR/README.md
#
# tfarrell@broadinstitute.org
# 20180404
#

library(reshape)
library(optparse)
suppressWarnings(suppressMessages(library(HaplotypR)))
suppressWarnings(suppressMessages(library(Biostrings)))
suppressWarnings(suppressMessages(library(ShortRead)))

# parse cmd-line options 
option_list <- list( 
    make_option(c("-o", "--output_dir"), 
                help="Directory where to save the output."),
    make_option(c("-p", "--amplicons_file"), 
                help="File with fwd/rev primers and reference seqs listed by amplicon."),
    make_option(c("-s", "--samples_dir"), 
                help="Directory with demultiplexed sample files."),
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE, 
                help="Run verbosely.")
)
opt = parse_args(OptionParser(option_list=option_list))

# prelims
setwd(opt$output_dir)
amplicons_file = opt$amplicons_file 
demultiplexed_samples = data.frame(FileR1 = list.files(opt$samples_dir, pattern="_R1.fastq.gz", full.names=TRUE),
                                   FileR2 = list.files(opt$samples_dir, pattern="_R2.fastq.gz", full.names=TRUE))
demultiplexed_samples$SampleID = lapply(demultiplexed_samples$FileR1, function(x) substr(basename(as.character(x)), 0, 13))
demultiplexed_samples$SampleName = lapply(demultiplexed_samples$FileR1, function(x) substr(basename(as.character(x)), 0, 13))
demultiplexed_samples$BarcodePair = lapply(seq(1, length(demultiplexed_samples$FileR1) * 2, 2), 
                                           function(x) { paste0("fwd_", x, "_rev_", x + 1) })
if (opt$verbose)  {
    cat("\ndemultiplexed_samples:\n")
    print(head(demultiplexed_samples))
} 

# demultiplex amplicons 
amplicon_df = read.delim(amplicons_file)
demultiplexed_amplicons_dir = "demultiplexed_amplicons/"
if (!dir.exists(demultiplexed_amplicons_dir)) { 
    cat("\ndemultiplexing amplicons...")
    dir.create(demultiplexed_amplicons_dir)
    demultiplexed_amplicons = lapply(2:dim(amplicon_df)[1], function(i) {
        mid = as.character(amplicon_df[i, "MarkerID"])
        adapter_fwd = as.character(amplicon_df[i, "Forward"])
        adapter_rev = as.character(amplicon_df[i, "Reverse"])
        r = lapply(seq_along(demultiplexed_samples$FileR1), function(j) { 
            output_file = file.path(demultiplexed_amplicons_dir, 
                                    sub("R1\\.fastq.gz", mid,
                                        basename(as.character(demultiplexed_samples$FileR1)[j])))
            removePrimer(as.character(demultiplexed_samples$FileR1)[j], 
                         as.character(demultiplexed_samples$FileR2)[j],
                         output_file, adapter_fwd, adapter_rev, 
                         max.mismatch=2, with.indels=FALSE)
        })
        cbind(BarcodePair=as.character(demultiplexed_samples$BarcodePair),
              MarkerID=mid, do.call(rbind, r))
    })
    demultiplexed_amplicons = do.call(rbind, demultiplexed_amplicons)
    demultiplexed_amplicons = merge.data.frame(demultiplexed_samples[,c("SampleID","SampleName","BarcodePair")],
                                               demultiplexed_amplicons, by="BarcodePair")
} else { 
    cat("\nretrieving demultiplexed amplicons...")
    demultiplexed_amplicons = data.frame(FileR1=character(),
                                         FileR2=character(),
                                         MarkerID=character(),
                                         SampleID=character(), 
                                         SampleName=character(),
                                         BarcodePair=character())
    for (amplicon in amplicon_df$MarkerID) { 
        df = data.frame(MarkerID=amplicon,
                        FileR1=list.files(demultiplexed_amplicons_dir, 
                                          pattern=paste0(amplicon, "_F"), full.names=TRUE),
                        FileR2=list.files(demultiplexed_amplicons_dir, 
                                          pattern=paste0(amplicon, "_R"), full.names=TRUE))
        df$SampleID = lapply(df$FileR1, function(x) substr(basename(as.character(x)), 0, 13))
        df$SampleName = df$SampleID
        df$BarcodePair = demultiplexed_samples[demultiplexed_samples$SampleID %in% df$SampleID,"BarcodePair"]
        demultiplexed_amplicons = rbind(demultiplexed_amplicons, df)
    }
}
if (opt$verbose) { 
    cat("\ndemultiplexed_amplicons:\n")
    print(head(demultiplexed_amplicons))
} 

# trim and merge paired-end reads 
read_lens_fwd = list(csp=185, sera2=175)
read_lens_rev = list(csp=103, sera2=84)
source("/Users/tfarrell/Tools/HaplotypR/R/processReads.R")
postfix = list()
for (amplicon in amplicon_df$MarkerID) { 
    postfix[[amplicon]] <- sprintf("_bind%.0f_%.0f", read_lens_fwd[[amplicon]], read_lens_rev[[amplicon]])
} 
processed_reads_dir = file.path("processed_reads")
if (!dir.exists(processed_reads_dir)) { 
    cat("\nprocessing reads...\n")
    dir.create(processed_reads_dir)
    # merge paired-end reads
    processed_reads = bindAmpliconReads(as.character(demultiplexed_amplicons$FileR1), 
                                        as.character(demultiplexed_amplicons$FileR2), 
                                        processed_reads_dir, read1Length=read_lens_fwd, read2Length=read_lens_rev,
                                        marker=as.character(demultiplexed_amplicons$MarkerID))
    processed_reads = cbind(demultiplexed_amplicons[,c("SampleID","SampleName","BarcodePair","MarkerID")], 
                            processed_reads)
} else {
    cat("\nretrieving processed reads...")
    processed_reads = demultiplexed_amplicons[,c("SampleID","SampleName","BarcodePair","MarkerID")]
    processed_reads = processed_reads[order(as.character(processed_reads[,"SampleID"])),]
    rownames(processed_reads) = 1:length(rownames(processed_reads))
    processed_reads[["ReadFile"]] = list.files(processed_reads_dir, full.names=TRUE)
}
if (opt$verbose) { 
    cat("\nprocessed_reads:\n")
    print(head(processed_reads))
} 

# compute mismatch rate and call SNPs
source("/Users/tfarrell/Tools/HaplotypR/R/callGenotype.R")
source("/Users/tfarrell/Tools/HaplotypR/R/callHaplotype.R")
cat("\ncomputing mismatch rates and calling SNPs...\n")
min_mismatch_rate = 0.5
min_genotype_occurence = 2
ref_seq = DNAStringSet(as.character(amplicon_df$ReferenceSequence))
names(ref_seq) = amplicon_df$MarkerID
snps_file = "snps.tsv"
if (!file.exists(snps_file)) { 
    snp_list = list()
    snp_df = data.frame(Chr=character(), Pos=integer(), Ref=character(), Alt=character())
    for (marker in amplicon_df$MarkerID) {
        # compute mismatch
        seq_errs = calculateMismatchFrequencies(as.character(processed_reads[processed_reads$MarkerID == marker, "ReadFile"]), 
                                                ref_seq[marker], method="pairwiseAlignment", minCoverage=100L)
        names(seq_errs) = processed_reads[processed_reads$MarkerID == marker, "SampleID"][seq(5)]
        seq_err = do.call(cbind, lapply(seq_errs, function(l) { l[,"MisMatch"]/l[,"Coverage"] }))
        # call SNPs
        possible_snp = callGenotype(seq_err, minMismatchRate=min_mismatch_rate, 
                                    minReplicate=min_genotype_occurence)
        snp_ref = unlist(lapply(possible_snp, function(snp) { as.character(subseq(ref_seq[marker], 
                                                                                  start=snp, width=1)) }))
        snps = cbind(Chr=marker, Pos=possible_snp, Ref=snp_ref, Alt="N")
        snp_list[[marker]] = snps
        if (length(possible_snp) > 0) {
            snp_df = rbind(snp_df, snps)
        }
    }
    write.table(snp_df, snps_file, quote=FALSE)
} else { 
    snp_list = list()
    snp_df = read.table(snps_file) 
    for (marker in amplicon_df$MarkerID) { 
        snp_list[[marker]] = snp_df[snp_df$Chr == marker,]
    }
}
if (opt$verbose) { 
    cat("\nsnps_df:\n")
    print(head(snp_df))
    print(snp_list)
} 

# call haplotypes
cat("\ncalling haplotypes...\n")
detection_limit = 0.01
min_haplotype_coverage = 3
min_haplotype_occurence = 2
min_sample_coverage = 25 
# call final haplotypes
haplotypes = createFinalHaplotypeTable(outputDir=opt$output_dir, 
                                       sampleTable=processed_reads, 
                                       snpLst=snp_list, refSeq=ref_seq, 
                                       postfix=postfix, markerTab=amplicon_df, 
                                       minReplicate=min_haplotype_occurence,
                                       minHaplotypCoverage=min_haplotype_coverage,
                                       minSampleCoverage=min_sample_coverage,
                                       detectability=detection_limit, include_seq=TRUE, 
                                       verbose=TRUE, just_contingency_table=FALSE)
cat("\ndone.\n\n")