vmatchPattern2_MC = function(pattern,subject,max.mismatch=0L,with.indels = T,ncpu = ncpu){
    hits = mclapply(split(subject,cut(1:length(subject),ncpu)),function(x){
        vmatchPattern2(pattern,x,max.mismatch = max.mismatch,with.indels = with.indels)
},mc.cores = ncpu)
    hits = Reduce("c",hits)
    return(hits)
}

vmatchPattern2 <- function(pattern, subject,
                              max.mismatch=0, min.mismatch=0,
                              with.indels=FALSE, fixed=TRUE,
                              algorithm="auto")
{
     if (!is(subject, "XStringSet"))
         subject <- Biostrings:::XStringSet(NULL, subject)
     algo <- Biostrings:::normargAlgorithm(algorithm)
     if (Biostrings:::isCharacterAlgo(algo))
         stop("'subject' must be a single (non-empty) string ",
              "for this algorithm")
     pattern <- Biostrings:::normargPattern(pattern, subject)
     max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
     min.mismatch <- Biostrings:::normargMinMismatch(min.mismatch,
max.mismatch)
     with.indels <- Biostrings:::normargWithIndels(with.indels)
     fixed <- Biostrings:::normargFixed(fixed, subject)
     algo <- Biostrings:::selectAlgo(algo, pattern,
                                     max.mismatch, min.mismatch,
                                     with.indels, fixed)
     C_ans <- .Call2("XStringSet_vmatch_pattern", pattern, subject,
                     max.mismatch, min.mismatch,
                     with.indels, fixed, algo,
                     "MATCHES_AS_RANGES",
                     PACKAGE="Biostrings")
     unlisted_ans <- IRanges(start=unlist(C_ans[[1L]],
use.names=FALSE),
                             width=unlist(C_ans[[2L]],
use.names=FALSE))
     relist(unlisted_ans, C_ans[[1L]])
}

scan_adapter = function(fqs,adapter,mismatch,ncpu){
    adapter.loc.df = mclapply(1:length(fqs),function(i){
        fq = as.character(fqs[[i]])
        fqname = sapply(str_split(names(fqs[i])," "),`[`,1)
        res = edlibR::align(adapter,fq,mode = "HW",task="path",k = mismatch)
        loc = res$locations
        if (length(loc) == 0){
            return(data.frame(fqname = character(),start = numeric(),end = numeric()))
        }
        start = min(sapply(loc,`[`,1)) + 1
        end = max(sapply(loc,`[`,2)) + 1
        adapter.loc = data.frame(fqname = fqname,start = start,end = end)
        return(adapter.loc)
    },mc.cores = ncpu)
    adapter.loc.df2 = do.call("bind_rows",adapter.loc.df)
    return(adapter.loc.df2)
}

plot_adapter_num = function(df,title){
    set.seed(43)
    palette <- distinctColorPalette(450)
    adapter.list.edlibR.mismatch.count = df %>% group_by(fqname) %>% summarise(adapter_count = n())
    df2 = df %>% group_by(fqname) %>%
        summarise(adapter_chain = paste(type,collapse = ","))
    df2$classification = unlist(lapply(df2$adapter_chain,classify_read))
    adapter.list.edlibR.mismatch.count.table = data.frame(table(adapter.list.edlibR.mismatch.count$adapter_count,df2$classification))
    colnames(adapter.list.edlibR.mismatch.count.table) = c("adapter_count","adapter_classification","count")
    adapter.list.edlibR.mismatch.count.table = adapter.list.edlibR.mismatch.count.table %>% arrange(adapter_classification)
    names(palette) = all.classification
    ggplot(adapter.list.edlibR.mismatch.count.table, aes(fill=adapter_classification, y=count, x=adapter_count)) + 
        geom_bar(position="stack", stat="identity") + scale_fill_manual(values = palette) +
    theme_bw() + guides(fill=guide_legend(ncol=2)) + ggtitle(title)
}

adapter_num = function(df,ncpus){
    adapter.list.edlibR.mismatch.count = df %>% group_by(fqname) %>% summarise(adapter_count = n())
    df2 = df %>% group_by(fqname) %>%
        summarise(adapter_chain = paste(type,collapse = ","))
    df2$classification = unlist(mclapply(df2$adapter_chain,classify_read,mc.cores = ncpus))
    return(df2)
}

order_suffix = function(x){
    paste(sort(str_split(x,pattern = "-")[[1]]),collapse  = "-")
}

classify_read = function(adapter_chain,adapters_graph = adapter_graph){
    adapters = str_split(adapter_chain,pattern = ",")[[1]]
    all.forward = all(!grepl(adapters,pattern = "*rc$"))
    all.reverse = all(grepl(adapters,pattern = "*rc$"))
    #This no.rc means that the read has no palindrome sequence.
    rc_mix = !(all.forward | all.reverse)
    forward.index = which(!grepl(adapters,pattern = "*rc$"))
    reverse.index = grep(adapters,pattern = "*rc$")

    rc_mix = ifelse(rc_mix,ifelse(all(diff(forward.index)==1) & all(diff(reverse.index)==1),"RC-successive","RC-mosaic"),"No-RC")
    forward.adapters = adapters[forward.index]
    reverse.adapters = adapters[reverse.index]
    
    if((length(forward.adapters)+length(reverse.adapters))>1){
        if(length(forward.adapters) == 1 & length(reverse.adapters) == 1){
            suffix = "Singleton-Singleton"
        }else if(length(forward.adapters) <= 1){
            reverse.adapters = data.frame(from = reverse.adapters[1:(length(reverse.adapters)-1)],
                                      to = reverse.adapters[2:length(reverse.adapters)])
            reverse.dist = apply(reverse.adapters,1,function(reverse.adapter){
                test_paths(adapters_graph,reverse.adapter[1],reverse.adapter[2])
            })
            reverse.adapter_sequence = ifelse(any(is.infinite(reverse.dist)),"Wrong direction",ifelse(all(reverse.dist==1),"Successive","Jump"))
            forward.adapter_sequence = ifelse(length(forward.adapters)==0,"NA","Singleton")
            suffix = paste(forward.adapter_sequence,reverse.adapter_sequence,sep = "-")
            suffix = order_suffix(suffix)
        }else if(length(reverse.adapters) <= 1){
            forward.adapters = data.frame(from = forward.adapters[1:(length(forward.adapters)-1)],
                                  to = forward.adapters[2:length(forward.adapters)])
            forward.dist = apply(forward.adapters,1,function(forward.adapter){
                test_paths(adapters_graph,forward.adapter[1],forward.adapter[2])
            })
            reverse.adapter_sequence = ifelse(length(reverse.adapters)==0,"NA","Singleton")
            forward.adapter_sequence = ifelse(any(is.infinite(forward.dist)),"Wrong direction",ifelse(all(forward.dist==1),"Successive","Jump"))
            suffix = paste(forward.adapter_sequence,reverse.adapter_sequence,sep = "-")
            suffix = order_suffix(suffix)
        }else{
            forward.adapters = data.frame(from = forward.adapters[1:(length(forward.adapters)-1)],
                                  to = forward.adapters[2:length(forward.adapters)])
            forward.dist = apply(forward.adapters,1,function(forward.adapter){
                test_paths(adapters_graph,forward.adapter[1],forward.adapter[2])
            })

            reverse.adapters = data.frame(from = reverse.adapters[1:(length(reverse.adapters)-1)],
                                      to = reverse.adapters[2:length(reverse.adapters)])
            reverse.dist = apply(reverse.adapters,1,function(reverse.adapter){
                test_paths(adapters_graph,reverse.adapter[1],reverse.adapter[2])
            })
            forward.adapter_sequence = ifelse(any(is.infinite(forward.dist)),"Wrong direction",ifelse(all(forward.dist==1),"Successive","Jump"))
            reverse.adapter_sequence = ifelse(any(is.infinite(reverse.dist)),"Wrong direction",ifelse(all(reverse.dist==1),"Successive","Jump"))
            suffix = paste(forward.adapter_sequence,reverse.adapter_sequence,sep = "-")
            suffix = order_suffix(suffix)
        }
        final = paste(rc_mix,suffix,sep = ";")
    }else{
        final = "Singleton"
    }
}

test_paths <- function(g, from, to){
    ifelse(is.finite(c(shortest.paths(g, from,to,mode = "out"))),c(shortest.paths(g, from,  to,mode = "out")),Inf)
}

test_node_onpath = function(g,from,to,test_node){
    is.finite(test_paths(g,from,test_node)) & is.finite(test_paths(g,test_node,to))
}

library(randomcoloR)
set.seed(42)
palette <- distinctColorPalette(32)
adapters = c(paste("adapter",LETTERS[1:16],sep  = ""),paste("adapter",LETTERS[1:16],".rc",sep  = ""))
mismatch.plot = function(mismatch.df,title){
    adapter.df.mismatch.df2 = mismatch.df %>% group_by(fqname) %>% 
    summarise(adapter_pair = paste(type[1:(length(type) - 1)], type[2:length(type)], sep = "-")) %>%
    filter(!grepl(pattern = "NA",adapter_pair))
    adapter.df.mismatch.df2.table = table(adapter.df.mismatch.df2$adapter_pair)
    adapter.pairs = data.frame(pair = names(adapter.df.mismatch.df2.table),count = as.vector(adapter.df.mismatch.df2.table))
    adapter.pairs.split = str_split(adapter.pairs$pair,"-")
    from = sapply(adapter.pairs.split,`[`,1)
    to = sapply(adapter.pairs.split,`[`,2)
    count = as.vector(adapter.pairs$count)
    adapter.df.mismatch.df3 = data.frame(from_node = from,to_node = to,count = count)
    p = adapter.df.mismatch.df3 %>% 
      uncount(count) %>%
      make_long(from_node, to_node) %>%
      mutate(node = fct_relevel(node, rev(adapters)),
             next_node = fct_relevel(next_node, rev(adapters))) %>%
      ggplot(aes(x = x, 
                 next_x = next_x, 
                 node = node, 
                 next_node = next_node,
                 fill = factor(node))) +
      geom_alluvial()+scale_fill_manual(values = palette)+theme_void()+ggtitle(title)
}

test_two_adapter = function(x,y){
    if(grepl(x,pattern = "start")|grepl(y,pattern = "end")){
	classification = "Boundary"
	return(classification)
    }
    x.rc = grepl(x,pattern = "rc")
    y.rc = grepl(y,pattern = "rc")
    if(x.rc == y.rc){
        dist = test_paths(adapter_graph,from = x,to = y)
        if (is.infinite(dist)){
            classification = "Wrong-samedirection"
        }else if(dist == 0){
            classification = "Identical"
        }else if(dist == 1){
            classification = "Correct"
        }else{
            classification = "Jump"
        }
    }else{
        x = gsub(pattern = ".rc",replacement = "",x)
        y = gsub(pattern = ".rc",replacement = "",y)
        dist = min(test_paths(adapter_graph,from = x,to = y),test_paths(adapter_graph,from = y,to = x))
        if(dist == 0){
            classification = "RC-Identical"
        }else{
            classification = "RC-Wrong"
        }
    }
    return(classification)
}

get_fq_basename = function(fq){
    suffixes_to_remove <- c(".fq", ".fq.gz", ".fastq", ".fastq.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")
    suffixes_pattern <- paste(suffixes_to_remove, collapse = "|")
    suffixes_pattern <- paste0("(", suffixes_pattern, ")$")
    fq_basename <- sub(suffixes_pattern, "", fq)
}

subseq_MC = function(dna,starts,ends,ncpu){
    #Note that length of dna,starts,ends must be the same.
    seq = mclapply(split(1:length(dna),cut(1:length(dna),ncpu)),function(x){
        subseq(dna[x],starts[x],ends[x])
    },mc.cores = ncpu)
    seqs = Reduce("c",seq)
    return(seqs)
}

get_fqa_type = function(fqa){
    suffixes_fq = c(".fq$", ".fq.gz$", ".fastq$", ".fastq.gz$")
    suffixes_fa = c(".fa$", ".fa.gz$", ".fasta$", ".fasta.gz$")
    suffixes_fq <- paste(suffixes_fq, collapse = "|")
    suffixes_fa <- paste(suffixes_fa, collapse = "|")
    fileType = ifelse(grepl(pattern = suffixes_fq,fqa),"fastq",ifelse(grepl(pattern = suffixes_fa,fqa),"fasta"))
    return(fileType)
}
