library(stringr)

args = commandArgs(trailingOnly = TRUE)

CIDtsv = args[1] # the directory of the adapter location
bed = str_replace(CIDtsv,".tsv",".bed")

qname_CID.rename = read.table(CIDtsv,sep = ",",header = T)
qname_CID.rename2 = data.frame(name = qname_CID.rename$raw_name,start = qname_CID.rename$start,end = qname_CID.rename$end,name2 = qname_CID.rename$raw_name,score = "0",strand = qname_CID.rename$strand,block_start = qname_CID.rename$start,block_end = qname_CID.rename$end,item_RGB = "255,0,0")
write.table(qname_CID.rename2,bed,sep = "\t",quote = F,row.names = F,col.names = F)
