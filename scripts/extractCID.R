library(Biostrings)
library(dplyr)
library(parallel)
library(stringr)
library(qs)

args = commandArgs(trailingOnly = TRUE)
adapter_dir = args[1] # the directory of the adapter location
max_mismatch = args[2]
ncpus = as.numeric(args[3])
fq = args[4]
CIDadapterids = args[5]
leftAdapter = args[6] # related to RNA(transcript) strand
rightAdapter = args[7] # related to RNA(transcript) strand

#extract fasta names
fqnames_file = paste0(fq,".fqnames")
cmd = paste0("zcat ",fq," |grep '^>'|sed 's/>//g' > ",fqnames_file)
system(cmd)

fqnames = read.table(fqnames_file,header = F) %>% pull(V1)

CIDadapterids = read.table(CIDadapterids,header = FALSE) %>% pull(V1)

mismatch_pattern = paste0(".mismatch",max_mismatch)

CIDadapter.files = list.files(adapter_dir,
                             pattern = paste0(".*splitted.*.*",mismatch_pattern,".*qs"),
                            full.names = TRUE)

adapter.list.mismatch = mclapply(CIDadapter.files,function(file){
    df = as.data.frame(qread(file))
    type2 = str_extract(file,pattern = paste(paste0(CIDadapterids,".*"),collapse="|"))
    type2 = str_replace(str_replace(type2,mismatch_pattern,""),".qs","")
    df$type = type2
    df = df %>% select(-group_name) %>% mutate(fqname = fqnames[group]) %>% select(-group)
    return(df)
},mc.cores=100)

adaptertypes = str_replace(str_replace(str_extract(CIDadapter.files,pattern = paste(paste0(CIDadapterids,".*"),collapse="|")),mismatch_pattern,""),".qs","")

names(adapter.list.mismatch) = adaptertypes

adapter_left_norc = adapter.list.mismatch[[leftAdapter]]
adapter_right_norc = adapter.list.mismatch[[rightAdapter]]
adapter_norc = full_join(adapter_left_norc,adapter_right_norc,by = "fqname") %>% select(fqname,end.x,start.y)
colnames(adapter_norc) = c("fqname","start","end")
adapter_norc$start = adapter_norc$start
adapter_norc$end = adapter_norc$end - 1
adapter_norc$type = ifelse(is.na(adapter_norc$start),"rightAdapter",ifelse(is.na(adapter_norc$end),"leftAdapter","both"))
adapter_norc$start[which(is.na(adapter_norc$start))] = pmax(adapter_norc$end[which(is.na(adapter_norc$start))] - 25,0)
adapter_norc$end[which(is.na(adapter_norc$end))] = adapter_norc$start[which(is.na(adapter_norc$end))] + 25
adapter_norc = adapter_norc %>% mutate(width = end - start) %>% filter(width>=20 & width<=30)

adapter_left_rc = adapter.list.mismatch[[paste0(leftAdapter,".rc")]]
adapter_right_rc = adapter.list.mismatch[[paste0(rightAdapter,".rc")]]
adapter_rc = full_join(adapter_left_rc,adapter_right_rc,by = "fqname") %>% select(fqname,start.x,end.y)
colnames(adapter_rc) = c("fqname","end","start")
adapter_rc$end = adapter_rc$end - 1 
adapter_rc$start = adapter_rc$start
adapter_rc$type = ifelse(is.na(adapter_rc$start),"leftAdapterRC",ifelse(is.na(adapter_rc$end),"rightAdapterRC","bothRC"))
adapter_rc$start[which(is.na(adapter_rc$start))] = pmax(adapter_rc$end[which(is.na(adapter_rc$start))] - 25,0)
adapter_rc$end[which(is.na(adapter_rc$end))] = adapter_rc$start[which(is.na(adapter_rc$end))] + 25
adapter_rc = adapter_rc %>% mutate(width = end - start) %>% filter(width>=20 & width<=30)

CID_loc = bind_rows(adapter_norc,adapter_rc)

out_file = paste0(adapter_dir,"/CID_loc.tsv")

write.table(CID_loc,file = out_file,col.names = T,row.names = F,quote = F)



