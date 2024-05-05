suppressMessages({
  library(Biostrings)
  library(ShortRead)
  library(stringr)
  library(dplyr)
  library(qs)
  library(GenomicRanges)
  library(parallel)
  source("scripts/utils.R")
})

args = commandArgs(trailingOnly = TRUE)
fq = args[1]
adapter_id = args[2]
adapter_seq = args[3]
maxmismatch = as.numeric(args[4])
ncpu = as.numeric(args[5])
outdir = args[6]

sample_base = get_fq_basename(basename(fq))

sample_type = get_fqa_type(basename(fq))

message(sample_type)

adapter_loc.all = IRangesList()
adapter_loc.rc.all = IRangesList()
width.df = data.frame(read = character(),width = numeric())

message(paste0("Processing adapter:",adapter_id))
message(paste0("Processing fastq:",fq))

adapter.file = paste0(outdir,"/",sample_base,"_",adapter_id,".mismatch",maxmismatch,".qs")
message(adapter.file)
adapter.rc.file = paste0(outdir,"/",sample_base,"_",adapter_id,".mismatch",maxmismatch,".rc.qs")
message(adapter.rc.file)
width.file = paste0(outdir,"/",sample_base,"_width",".qs")
message(width.file)

fq <- open_input_files(fq)

count = 0
while (TRUE) {
    ## Load 3e5 records at a time. Each new call to readDNAStringSet()
    ## picks up where the previous call left.
    chunk = 3e5
    count = count + chunk
    dna <- readDNAStringSet(fq, nrec = chunk,format = sample_type)
    if (length(dna) == 0L)
        break
    cat("processing reads:", count, "...\n")

    adapter_loc = vmatchPattern2_MC(adapter_seq,dna,
      max.mismatch = maxmismatch,with.indels = T,ncpu = ncpu)
    adapter_loc.rc = vmatchPattern2_MC(reverseComplement(DNAString(adapter_seq)),
      dna,max.mismatch = maxmismatch,with.indels = T,ncpu = ncpu)

    adapter_loc.all = c(adapter_loc.all,adapter_loc)
    adapter_loc.rc.all = c(adapter_loc.rc.all,adapter_loc.rc)

    width.df = rbind(width.df,data.frame(read = names(dna),width = width(dna)))
}

width.df$read = unlist(mclapply(strsplit(width.df$read," "), `[`, 1,mc.cores = ncpu))

qsave(adapter_loc.all,file = adapter.file)
qsave(adapter_loc.rc.all,file = adapter.rc.file)
qsave(width.df,file = width.file)
