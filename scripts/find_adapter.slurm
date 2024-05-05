#! /bin/bash

#SBATCH -J find_adapter
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 100
#SBATCH --array=0-15

source activate jupyter

fastq=$1
mismatch=$2
outdir=$3
adapterids=$4
adapterseqs=$5

readarray -t adapterids < ${adapterids}

readarray -t adapterseqs < ${adapterseqs}

adapterid=${adapterids[$SLURM_ARRAY_TASK_ID]}
adapterseq=${adapterseqs[$SLURM_ARRAY_TASK_ID]}

Rscript scripts/find_adapter.R \
 ${fastq} ${adapterid} ${adapterseq} ${mismatch} 100 ${outdir}
