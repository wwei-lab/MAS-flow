#! /bin/bash

#SBATCH -J segment_classfication 
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 100

source activate jupyter

outdir=$1
fastq=$2
mismatch=$3
width_file=`ls ${outdir}/*width.qs`

echo ${width_file}

Rscript scripts/segment-classify.R ${outdir} ${width_file} ${mismatch} 100


