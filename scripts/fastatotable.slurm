#! /bin/bash

#SBATCH -J fatotable
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 2

cat `cat CIDfasta.txt`|awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print $1,$3}' > all_CID.fasta.tsv 
