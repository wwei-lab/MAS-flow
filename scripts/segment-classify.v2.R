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

tsv.Jump = tsv %>% filter(grepl(fqname,pattern = "Jump"))

tsv.no_split = tsv %>% filter(grepl(fqname,pattern = "no-split"))

##############Deal with Jump-adapter reads:01-filter##############
adapterPairs = mclapply(str_split(tsv.Jump$fqname,";"),function(x){x[length(x)-1]},mc.cores = ncpus)

adapterPairs = unlist(str_extract_all(adapterPairs,pattern = "adapter."))
start.adapter = adapterPairs[seq(1,length(adapterPairs),2)]
end.adapter = adapterPairs[seq(2,length(adapterPairs),2)]
tsv.Jump$start_adapter = start.adapter
tsv.Jump$end_adapter = end.adapter

segment_width = mclapply(str_split(tsv.Jump$fqname,";"),function(x){x[length(x)-2]},mc.cores = ncpus)
segment_width = mclapply(str_split(segment_width,":"),function(x){x[3]},mc.cores = ncpus)
st_ed = unlist(str_split(segment_width,"-"))
start_pos = st_ed[seq(1,length(st_ed),2)]
end_pos = st_ed[seq(2,length(st_ed),2)]
segment_width = as.numeric(end_pos) - as.numeric(start_pos) + 1

tsv.Jump$segment_width = segment_width

keep = which(unlist(mclapply(1:nrow(tsv.Jump),function(x){
    test_node_onpath(g = adapter_graph,tsv.Jump$start_adapter[x],tsv.Jump$end_adapter[x],test_node = tsv.Jump$type[x])
},mc.cores = ncpus)))

tsv.Jump.keep = tsv.Jump[keep,]

#Deal with Jump-adapter reads:02-reformt it to trim-ready data frame
tsv.Jump.keep.middle = data.frame(fqname = tsv.Jump.keep$fqname,start = tsv.Jump.keep$start,
                        end = tsv.Jump.keep$end,type = tsv.Jump.keep$type)

#magic 1-bp adapter
tsv.Jump.keep.boundary.startadapter = unique(data.frame(fqname = tsv.Jump.keep$fqname,
                        start = 0,end = 0,type = tsv.Jump.keep$start_adapter))


tsv.Jump.keep.boundary.endadapter = unique(data.frame(fqname = tsv.Jump.keep$fqname,
                        start = tsv.Jump.keep$segment_width+1,end = tsv.Jump.keep$segment_width+1,type = tsv.Jump.keep$end_adapter))

tsv.Jump.keep.boundary = rbind(tsv.Jump.keep.boundary.startadapter,tsv.Jump.keep.boundary.endadapter)

tsv.Jump.keep = rbind(tsv.Jump.keep.boundary,tsv.Jump.keep.middle)

##############Deal with Jump-adapter done##############

##############Deal with no-split-adapter###############
tsv.no_split.keep = tsv.no_split %>% select(fqname,start,end,type)

segment_width = mclapply(str_split(tsv.no_split.keep$fqname,";"),function(x){x[length(x)-2]},mc.cores = ncpus)
segment_width = mclapply(str_split(segment_width,":"),function(x){x[3]},mc.cores = ncpus)
st_ed = unlist(str_split(segment_width,"-"))
start_pos = st_ed[seq(1,length(st_ed),2)]
end_pos = st_ed[seq(2,length(st_ed),2)]
segment_width = as.numeric(end_pos) - as.numeric(start_pos) + 1

tsv.no_split.keep$width = segment_width

width.df = tsv.no_split.keep %>% select(fqname,width) %>% unique()

end.df = data.frame(start = width.df$width+1,end = width.df$width+1,type = "end",fqname = width.df$fqname,width = 1) %>%
select(fqname,start,end,type,width)

start.df = data.frame(start = 0,end = 0,type = "start",fqname = width.df$fqname,width = 1)%>%
select(fqname,start,end,type,width)

tsv.no_split.keep = rbind(rbind(tsv.no_split.keep,start.df),end.df) %>% arrange(fqname,start)
##############Deal with no-split-adapter done############

tsv2 = bind_rows(tsv.Jump.keep,tsv.no_split.keep)

tsv2 = tsv2 %>% arrange(fqname,start) %>% group_by(fqname)

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
       group_by(fqname) %>% mutate(segment = paste(fqname,row_number(),sep = "|"))

mismatch_out2 = paste0(adapter_dir,"/adapter.mismatch",max_mismatch,".moresplit.tosplit.tsv")

write.table(tsv3,mismatch_out2,quote = F,row.names = F,col.names = F)

message("To-split dataframe done!")






