import pysam
import pandas as pd
import numpy as np
import sys

inbam = sys.argv[1]
outdir = sys.argv[2]
outbam = inbam.replace(".bam",".rename.bam")

infile = pysam.AlignmentFile(inbam, "r")
outfile = pysam.AlignmentFile(outbam, "w", template=infile)
count = 0
id_CID_out = outdir + "/qname_CID.tsv"
aln = id_CID_out
aln_out = outdir + "/qname_CID.rename.tsv"
CID_loc_file = outdir + "/CID_loc.tsv"

id_CID_out_file = open(id_CID_out,"w")

for s in infile:
    if count%100000 == 0:
        print("processing:",count)

    qstart = s.query_alignment_start

    qend = s.query_alignment_end
    oldname = s.query_name
    s.query_name = s.query_name + "|" + "ts-te:" + str(qstart) + "-" + str(qend)
    
    outfile.write(s)
    outline = oldname+"\t"+s.query_name+"\t"+str(qstart)+"\t"+str(qend)+"\n"
    id_CID_out_file.write(outline)
    count += 1

id_CID_out_file.close()
infile.close()
outfile.close()

def if_substring(a,pattern):
    b = np.array([True if i.find(pattern)!= -1 else False for i in a])
    return(b)

CID_loc = pd.read_table(CID_loc_file,sep=" ")
aln_pd  = pd.read_table(aln,sep="\t",names = ["raw_name","new_name","qstart","qend"])
aln_pd2 = aln_pd.merge(CID_loc,how = "left",left_on = "raw_name",right_on = "fqname")
types = aln_pd2["type"].to_numpy().astype("str")
qstart = aln_pd2["qstart"].to_numpy()
qend = aln_pd2["qend"].to_numpy()
start = aln_pd2["start"].to_numpy().astype("int") # for 0-base which can be used to generate bed
end = aln_pd2["end"].to_numpy().astype("int")
aln_pd2["strand"] =  np.where(if_substring(types,"RC"),"-","+")
aln_pd2["distance"] =  np.where(if_substring(types,"RC"),qstart - end,start - qend)
aln_pd2["distance2"] = abs(aln_pd2["distance"] - 50)
aln_pd2 = aln_pd2[aln_pd2['distance2'] == aln_pd2.groupby('new_name')['distance2'].transform('min')]
aln_pd2 = aln_pd2[aln_pd2["distance2"]<=80]
ranges = (aln_pd2["start"]+1).astype(int).astype(str) + "-" + aln_pd2["end"].astype(int).astype(str)
aln_pd2["CIDname"] = aln_pd2["raw_name"] + "_" + ranges + ":" + aln_pd2["strand"]

aln_pd2.to_csv(aln_out, index=False)
