#! /bin/bash

#SBATCH -J find_CIDoligo
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 100
#SBATCH --array=0-2


source activate jupyter

readarray -t adapterids < data/adapterids.head.txt

readarray -t adapterseqs < data/adapterseqs.head.txt

adapterid=${adapterids[$SLURM_ARRAY_TASK_ID]}
adapterseq=${adapterseqs[$SLURM_ARRAY_TASK_ID]}

fastq=$1
mismatch=$2
outdir=$3

Rscript scripts/find_adapter.R \
 ${fastq} ${adapterid} ${adapterseq} ${mismatch} 100 ${outdir}
