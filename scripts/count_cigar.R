library(Rsamtools)
library(qs)
library(parallel)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
bam = args[1]
outdir = args[2]
Q = args[3]

readBam <- function(bamfile,chunk,Q) {

  bf <- open(BamFile(bamfile, yieldSize=chunk))
  para = ScanBamParam(tag = c("NM"),what = c("qname","flag", "rname", 
                                             "pos", "mapq","cigar"),mapqFilter = Q)
  ops <- GenomicAlignments::CIGAR_OPS

  count = 0
  bam.df.all = data.frame()
  while (TRUE){
    bam = scanBam(bf,param = para)[[1]]
    if(length(bam[[1]]) == 0)
      break
    count = count + chunk
    cat("processing:", count, "\n")
    widths <- GenomicAlignments::explodeCigarOpLengths(bam[["cigar"]], ops = ops)
    keep.ops <- GenomicAlignments::explodeCigarOps(bam[["cigar"]], ops = ops)
    explodedcigars <- IRanges::CharacterList(relist(paste0(unlist(widths), 
                                                         unlist(keep.ops)), widths))
    bam.df = data.frame(qname = bam[["qname"]],flag = bam[["flag"]],
              rname = bam[["rname"]],pos = bam[["pos"]],mapq = bam[["mapq"]],NM = bam[["tag"]][["NM"]])
    for (opt in setdiff(ops, "=")) {
      bam.df[,paste0("nbr",opt)] = 
        unlist(mclapply(explodedcigars, function(cg) sum(as.numeric(gsub(paste0(opt, "$"), "", cg)),na.rm = TRUE),mc.cores = 100))
    }
    bam.df.all = rbind(bam.df.all,bam.df)
  }
  return(bam.df.all)
}

tsv = readBam(bam,3e6,Q)

tsv = tsv %>% mutate(base_accuracy = (nbrM+nbrI+nbrD-NM)/(nbrM+nbrI+nbrD))

outfile = paste0(outdir,"/","bamcount.mapq",Q,".tsv")

qsave(tsv,file = outfile)
