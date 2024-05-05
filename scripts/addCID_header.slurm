#! /bin/bash

#SBATCH -J getCIDtable
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 1
#SBATCH -c 1
#SBATCH --array=0-8

readarray -t outdirs < CIDtbltobed.tsv

outdir=${outdirs[$SLURM_ARRAY_TASK_ID]}

bam=`ls ${outdir}/*.bam`

python scripts/addCID_header.py ${bam} ${outdir}
