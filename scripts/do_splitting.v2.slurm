#! /bin/bash

#SBATCH -J do_splitting
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 100

source activate jupyter

outdir=$1
fa=$2
mismatch=$3

Rscript scripts/do_splitting.v2.R ${outdir} ${mismatch} 100 ${fa}
