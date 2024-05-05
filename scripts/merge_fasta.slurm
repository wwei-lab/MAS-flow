#! /bin/bash

#SBATCH -J merge_fasta
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 1
#SBATCH -c 1

moresplit_trim_fa2=$1
nomoresplit_trim_fa=$2
trim_fa=$3

cat ${moresplit_trim_fa2} ${nomoresplit_trim_fa} > `dirname ${trim_fa}`/`basename ${trim_fa} .trimmed.fasta.gz`.splitted.fasta.gz
