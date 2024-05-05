library(dplyr)
library(parallel)
library(stringr)
library(Biostrings)
library(qs)
source("scripts/utils.R")
source("scripts/adapter_graph.R")

##################mismatch dataframe done####################################
args = commandArgs(trailingOnly = TRUE)
adapter_dir = args[1] # the directory of the adapter location
width_file = args[2] # width file
max_mismatch = args[3]
ncpus = args[4]
mismatch_out = paste0(adapter_dir,"/adapter.list.mismatch",max_mismatch,".df.tsv")

width.df = qread(width_file)

fqnames = width.df$read

adapter.files.mismatch = list.files(adapter_dir,pattern = paste0("*adapter*.*mismatch",max_mismatch,".*qs$"),full.names = T)

adapter.list.mismatch = mclapply(adapter.files.mismatch,function(file){
    df = as.data.frame(qread(file))
    mismatchpattern = paste0(".","mismatch",max_mismatch) 
    type2 = str_replace(str_replace(str_extract(file,
                        pattern = "adapter.*"),mismatchpattern,""),".qs","")
    df$type = type2
    df = df %>% select(-group_name) %>% mutate(fqname = fqnames[group]) %>% select(-group)
    return(df)
},mc.cores = ncpus)

adapter.list.mismatch= do.call("bind_rows",adapter.list.mismatch)
adapter.list.mismatch.df = adapter.list.mismatch %>% arrange(fqname,start)%>%
filter(grepl(type,pattern = "adapter"))
write.table(adapter.list.mismatch.df,
            file = mismatch_out,quote = F,row.names = F,col.names =F)

message("Mismatch dataframe done!")

################mismatch dataframe done#################################

tsv = read.table(mismatch_out)
colnames(tsv) = c("start","end","width","type","fqname")

tsv2 = tsv %>% group_by(fqname) %>% filter(n() >= 2)

end.df =  width.df %>% filter(read %in% unique(tsv2$fqname))
end.df = data.frame(start = end.df$width+1,end = end.df$width+1,width = 1,type = "end",fqname = end.df$read)

start.df = data.frame(start = 0,end = 0,width = 1,type = "start",fqname = end.df$fqname)

tsv2 = rbind(rbind(tsv2,start.df),end.df) %>% arrange(fqname,start)

tsv2.head = tsv2 %>% filter(row_number() < n())
tsv2.tail = tsv2 %>% filter(row_number() > 1)

tsv2.head = tsv2.head %>%
    setNames(paste0('head.', names(.)))
tsv2.tail = tsv2.tail %>%
    setNames(paste0('tail.', names(.)))

segment_type = unlist(mclapply(1:nrow(tsv2.head),function(i){
    test_two_adapter(tsv2.head$head.type[i],
                     tsv2.tail$tail.type[i])
},mc.cores = ncpus))

tsv3 = data.frame(fqname = tsv2.head$head.fqname,
                  start = tsv2.head$head.end + 1,
                  end = tsv2.tail$tail.start - 1,
		  adapter = paste(tsv2.head$head.type,tsv2.tail$tail.type,sep = "-"),
                  segment_type = segment_type) %>% 
       group_by(fqname) %>% mutate(segment = paste(fqname,row_number(),sep = "-"))

mismatch_out2 = paste0(adapter_dir,"/adapter.mismatch",max_mismatch,".tosplit.tsv")

write.table(tsv3,mismatch_out2,quote = F,row.names = F,col.names = F)

message("To-split dataframe done!")
