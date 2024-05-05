library(Biostrings)
library(qs)
library(dplyr)
library(parallel)
library(stringr)
source("scripts/utils.R")

args = commandArgs(trailingOnly = TRUE)
adapter_dir = args[1] # the directory of the adapter location
max_mismatch = args[2]
ncpus = as.numeric(args[3])
fq = args[4]

message(adapter_dir)

sample_type = get_fqa_type(basename(fq))

hierachy_match = TRUE

overwrite = TRUE

fq_base = get_fq_basename(fq)

out_fa = paste0(fq_base,".trimmed.fasta.gz")

if(overwrite & file.exists(out_fa)){
    file.remove(out_fa)
}

width_file = list.files(adapter_dir,pattern = "*_width.qs",full.names = TRUE)

message(width_file)

width.df = qread(width_file)

mismatchpattern = paste0("mismatch",max_mismatch)
to_split_file = paste0(adapter_dir,"/adapter.",mismatchpattern,".tosplit.tsv")
to_split = read.table(to_split_file,header = F)
colnames(to_split) = c("read","start","end","adapterPairs","adapterType","segment")
to_split_fqnames = unique(to_split$read)

fq <- open_input_files(fq)

count = 0

while (TRUE) {
    ## Load 1e5 records at a time. Each new call to readDNAStringSet()
    ## picks up where the previous call left.
    chunk = 5e5
    count = count + chunk
    dna <- readDNAStringSet(fq, nrec = chunk,format = sample_type)
    if (length(dna) == 0L)
        break
    cat("processing reads:", count, "...\n")

    names(dna) = unlist(mclapply(names(dna),function(x) {str_split(x," ")[[1]][1]},mc.cores = ncpus))

    nosplit_fqnames = setdiff(names(dna),to_split_fqnames)
    no_split_dna = dna[nosplit_fqnames]

    names(no_split_dna) = paste(names(no_split_dna),
        paste("st_ed","A",paste(1,width(no_split_dna),sep = "-"),sep = ":"),"AD:A:no-split","AP:NA:NA",sep = ";")

    split_fqnames = intersect(names(dna),to_split_fqnames)
    to_split_target = to_split %>% filter(read %in% split_fqnames) %>% 
        filter((adapterType != "Boundary" & end - start + 1 >= 1) | (adapterType == "Boundary" & end - start + 1 >= 50))
    dna_splited = subseq(dna[to_split_target$read],to_split_target$start,to_split_target$end)
    st_ed_tag = paste("st_ed","A",paste(to_split_target$start,to_split_target$end,sep = "-"),sep = ":")
    AD_tag = paste("AD","A",to_split_target$adapterPairs,sep = ":") # tag for adapter Pairs
    AP_tag = paste("AP","A",to_split_target$adapterType,sep = ":") # tag for adapter Type
    #We add fastq comments after space of the header name,this will be pretained when adding -y in minimap2 alignments
    names(dna_splited) = paste(to_split_target$segment,st_ed_tag,AD_tag,AP_tag,sep = ";")

    dna_all = c(no_split_dna,dna_splited)

    writeXStringSet(dna_all, out_fa, compress=TRUE, format="fasta",append = ifelse(file.exists(out_fa),TRUE,FALSE))
}
