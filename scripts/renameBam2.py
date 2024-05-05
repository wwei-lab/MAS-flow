import pysam
import pandas as pd
import numpy as np
import sys

inbam = sys.argv[1]
outdir = sys.argv[2]
inbam2 = inbam.replace(".bam",".rename.bam")

outbam = inbam.replace(".bam",".rename2.bam")

infile = pysam.AlignmentFile(inbam2, "r")
outfile = pysam.AlignmentFile(outbam, "w", template=infile)
count = 0
aln_out = outdir + "/qname_CID.rename.tsv"

name_dict = dict()

count2 = 0

count3 = 0

with open(aln_out,"r")  as aln_out:
    for line in aln_out:
        count2 += 1
        if count2 == 1:
            continue
        lines = line.strip().split(",")
        new_name = lines[1]
        CID_name = lines[12]
        name_dict[new_name] = CID_name

for s in infile:
    if count%100000 == 0:
        print("processing:",count,"\t")
        print("processing success:",round(count3/(count+1),2))
    count += 1
    try:
        s.query_name = name_dict[s.query_name]
        count3 += 1
    except KeyError as e:
        pass
    outfile.write(s)

infile.close()
outfile.close()
