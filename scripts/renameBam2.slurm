#! /bin/bash

#SBATCH -J renameBam
#SBATCH -p q1
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mem=100G
#SBATCH --ntasks 10
#SBATCH --array=0-8

readarray -t bams < allrawbam.list.txt

readarray -t outdirs < outdirlist2.txt 

bam=${bams[$SLURM_ARRAY_TASK_ID]}

outdir=${outdirs[$SLURM_ARRAY_TASK_ID]}

python scripts/renameBam2.py ${bam} ${outdir}
