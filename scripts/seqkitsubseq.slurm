#! /bin/bash

#SBATCH -J seqkit_subseq
#SBATCH -p q1
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mem=50G
#SBATCH --ntasks 1
#SBATCH -c 100
#SBATCH --array=0-8 

readarray -t fqs < splitted.falist.txt

readarray -t beds < bedlist.txt 

fq=${fqs[$SLURM_ARRAY_TASK_ID]}

bed=${beds[$SLURM_ARRAY_TASK_ID]}

source activate seqkit

seqkit subseq -j 100 --bed ${bed} ${fq} > `dirname ${fq}`/`basename ${fq} .fasta.gz`.CID.fa
