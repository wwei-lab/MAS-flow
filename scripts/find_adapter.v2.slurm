#! /bin/bash

#SBATCH -J find_adapter_oligo
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 100
#SBATCH --array=0-15

fastq=$1
mismatch=$2
outdir=$3
CIDadapterids=$4
CIDadapterseqs=$5

source activate jupyter

readarray -t adapterids < ${CIDadapterids}

readarray -t adapterseqs < ${CIDadapterseqs}

adapterid=${adapterids[$SLURM_ARRAY_TASK_ID]}
adapterseq=${adapterseqs[$SLURM_ARRAY_TASK_ID]}

Rscript scripts/find_adapter.R \
 ${fastq} ${adapterid} ${adapterseq} ${mismatch} 100 ${outdir}
