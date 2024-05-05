import os
import pandas as pd

files = os.listdir("aln_cid_match.split/")

files.sort()

files = ["aln_cid_match.split/"+f for f in files]

match = pd.read_table(files[0])

for file in files[1:]:
    print(file)
    match2 = pd.read_table(file)
    match['CID'] = match['CID'].fillna(match2['CID'])

match.to_csv("aln_cid_match.split/match_all.v2.tsv",sep = ",")
