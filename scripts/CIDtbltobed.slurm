#! /bin/bash

#SBATCH -J CIDtbltobed
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 2
#SBATCH --mem=100G
#SBATCH -c 2
#SBATCH --array=0-8

readarray -t CIDtbls < CIDtbltobed.tsv

CIDtbl=${CIDtbls[$SLURM_ARRAY_TASK_ID]}

source activate jupyter

Rscript scripts/CIDtbltobed.R ${CIDtbl}
