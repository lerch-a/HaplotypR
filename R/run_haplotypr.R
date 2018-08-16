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
                help="Directory to save the output."),
    make_option(c("-p", "--amplicons_file"), 
                help="File with fwd/rev primers, reference seqs, fwd/rev read lengths and max indel threshold listed by amplicon."),
    make_option(c("-s", "--samples_dir"), 
                help="Directory with demultiplexed sample files split into R1/2\n\t\t(i.e. sample_demultiplexed/sample_*_R[12].fastq.gz)."),
    make_option(c("-t", "--trim_reads"), action="store_true", default=FALSE,
                help="If passed will look for fwd/rev read lengths in amplicons_file and trim reads to those lengths."),
    make_option(c("--min_mismatch"), default=0.05, 
                help="Minimum rate of mismatch between haplotype and reference sequence (default=0.05)."),
    make_option(c("--min_genotype_occurrence"), default=2, 
                help="Minimum # of samples for a valid genotype to be called in (default=2)."),
    make_option(c("--detection_limit"), default=0.01, 
                help="Minimum frequency for detecting a haplotype (default=0.01)."),
    make_option(c("--min_haplotype_coverage"), default=3, 
                help="Minimum coverage for haplotype to be recognized as valid (default=3)."),
    make_option(c("--min_haplotype_occurrence"), default=2, 
                help="Minimum # of samples for a valid haplotype to be called in (default=2)."),
    make_option(c("--min_sample_coverage"), default=2, 
                help="Minimum coverage for a sample to be recognized as valid (default=2)."),
    make_option(c("--return_full_haplotypes"), action="store_true", default=FALSE, 
                help=paste0("Pass this flag to return full haplotype sequences in the output\n",
                            "\t\t(as opposed to just sequence of alleles that occur at identified snp positions).")),
    make_option(c("--test_on"), default=0, 
                help="Test pipeline on the first $test_on samples in the sample-demultiplexed directory."),
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Run verbosely.")
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
if (opt$test_on > 0) { 
    demultiplexed_samples = demultiplexed_samples[seq(opt$test_on),]    
}
if (opt$verbose)  {
    cat(paste0("\ndemultiplexed_samples: ", as.character(nrow(demultiplexed_samples)), '\n'))
    print(head(demultiplexed_samples))
} 

# demultiplex amplicons 
amplicon_df = read.delim(amplicons_file)
demultiplexed_amplicons_dir = "demultiplexed_amplicons/"
if (!dir.exists(demultiplexed_amplicons_dir)) { 
    cat("\ndemultiplexing amplicons...\n")
    dir.create(demultiplexed_amplicons_dir)
    demultiplexed_amplicons = lapply(1:dim(amplicon_df)[1], function(i) {
        mid = as.character(amplicon_df[i, "MarkerID"])
        cat(paste0(mid, '...\n'))
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
    demultiplexed_samples$BarcodePair = as.character(demultiplexed_samples$BarcodePair)
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
postfix = list()
#source("/Users/tfarrell/Tools/HaplotypR/R/processReads.R")
if (opt$trim_reads) { 
    read_lens_fwd = as.integer(amplicon_df$read_lens_fwd)
    names(read_lens_fwd) = amplicon_df$MarkerID
    read_lens_rev = as.integer(amplicon_df$read_lens_rev)
    names(read_lens_rev) = amplicon_df$MarkerID
    for (amplicon in amplicon_df$MarkerID) { 
        postfix[[amplicon]] <- sprintf("_bind%.0f_%.0f", read_lens_fwd[[amplicon]], read_lens_rev[[amplicon]])
    } 
} else { 
    read_lens_fwd = list()
    read_lens_rev = list()
    for (amplicon in amplicon_df$MarkerID) { 
        postfix[[amplicon]] = ""
        #read_lens_fwd[[amplicon]] = c(NULL)
        #read_lens_rev[[amplicon]] = c(NULL)
    }
}
if (opt$verbose) { 
    cat("\npostfix and fwd/rev read lens:\n")
    print(postfix)
    print(read_lens_fwd)
    print(read_lens_rev)
}
processed_reads_dir = file.path("processed_reads")
if (!dir.exists(processed_reads_dir)) { 
    cat("\nprocessing reads...\n")
    dir.create(processed_reads_dir)
    # merge paired-end reads
    processed_reads = bindAmpliconReads(as.character(demultiplexed_amplicons$FileR1), 
                                        as.character(demultiplexed_amplicons$FileR2), 
                                        processed_reads_dir, read1Length=read_lens_fwd, read2Length=read_lens_rev,
                                        marker=as.character(demultiplexed_amplicons$MarkerID),
                                        postfix=postfix)
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
#source("/Users/tfarrell/Tools/HaplotypR/R/callGenotype.R")
#source("/Users/tfarrell/Tools/HaplotypR/R/callHaplotype.R")
cat("\ncomputing mismatch rates and calling SNPs...\n")
min_mismatch_rate = opt$min_mismatch 
min_genotype_occurrence = opt$min_genotype_occurrence 
ref_seq = DNAStringSet(as.character(amplicon_df$ReferenceSequence))
names(ref_seq) = amplicon_df$MarkerID
snps_file = "snps.tsv"
if (!file.exists(snps_file)) { 
    snp_list = list()
    snp_df = data.frame(Chr=character(), Pos=integer(), Ref=character(), Alt=character())
    for (marker in amplicon_df$MarkerID) {
        mismatch_freqs_file = paste0(as.character(marker), ".mismatch_freqs.tsv")
        if (!file.exists(mismatch_freqs_file)) { 
            if (opt$verbose) cat("\ncomputing mismatch freqs...\n")
            # compute mismatch freqs
            seq_errs = calculateMismatchFrequencies(as.character(processed_reads[processed_reads$MarkerID == marker, "ReadFile"]), 
                                                    as.character(ref_seq[[marker]]), method="pairwiseAlignment", minCoverage=100L)
            names(seq_errs) = processed_reads[processed_reads$MarkerID == marker, "SampleID"]
            seq_err = do.call(cbind, lapply(seq_errs, function(l) { l[,"MisMatch"]/l[,"Coverage"] }))
            write.table(seq_err, mismatch_freqs_file)
        } else { 
            seq_err = read.table(mismatch_freqs_file)    
        }
        # call SNPs
        if (opt$verbose) cat("\ncalling SNPs...\n")
        possible_snp = callGenotype(seq_err, minMismatchRate=min_mismatch_rate, 
                                    minReplicate=min_genotype_occurrence)
        snp_ref = unlist(lapply(possible_snp, function(snp) { as.character(subseq(as.character(ref_seq[[marker]]), 
                                                                                  start=snp, width=1)) }))
        snps = as.data.frame(cbind(Chr=marker, Pos=possible_snp, Ref=snp_ref, Alt="N"))
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

# call final haplotypes
cat("\ncalling haplotypes...\n")
max_indel_thresholds = amplicon_df$max_indel_threshold
names(max_indel_thresholds) = amplicon_df$MarkerID
haplotypes = createFinalHaplotypeTable(outputDir=opt$output_dir, 
                                       sampleTable=processed_reads, 
                                       snpLst=snp_list, refSeq=ref_seq, 
                                       postfix=postfix, markerTab=amplicon_df, 
                                       minReplicate=opt$min_haplotype_occurrence,
                                       minHaplotypCoverage=opt$min_haplotype_coverage,
                                       minSampleCoverage=opt$min_sample_coverage,
                                       detectability=opt$detection_limit, include_seq=TRUE, 
                                       verbose=TRUE, just_contingency_table=FALSE, 
                                       max_indel_thresholds=max_indel_thresholds,
                                       return_full_haplotypes=opt$return_full_haplotypes)
cat("\ndone.\n\n")