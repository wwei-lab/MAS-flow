#! /bin/bash

#SBATCH -J extractCID
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 100
#SBATCH --mem=100G
#SBATCH -c 1

source activate jupyter

outdir=$1
mismatch=$2
ncpus=$3
fq=$4
adapterid=$5
leftadapter=$6
rightadapter=$7

Rscript scripts/extractCID.R ${outdir} ${mismatch} ${ncpus} \
 ${fq} ${adapterid} ${leftadapter} ${rightadapter}
