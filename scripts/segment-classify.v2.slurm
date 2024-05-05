#! /bin/bash

#SBATCH -J segment_classfication 
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 100

source activate jupyter

mismatch=$1
outdir=$2
widthfile=`ls ${outdir}/*moresplit*width.qs`

echo ${outdir} ${widthfile} ${mismatch}

Rscript scripts/segment-classify.v2.R ${outdir} ${widthfile} ${mismatch} 100
