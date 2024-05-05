library(igraph)
library(dplyr)
library(parallel)
library(stringr)
library(randomcoloR)
library(ggplot2)
library(Biostrings)
source("scripts/utils.R")

tsv = read.table("data/FAX23635_pass_12b693ce_c74c2dc4//adapter.list.T44.mismatch3.df.tsv")
colnames(tsv) = c("start","end","width","type","fqname")
T47_new_adapter = adapter_num(tsv,80)

write.table(T47_new_adapter,file = "data/FAX23635_pass_12b693ce_c74c2dc4/adapter.T44.mismatch3.classification.tsv",quote = F,row.names = F,col.names = F,sep = "\t")
