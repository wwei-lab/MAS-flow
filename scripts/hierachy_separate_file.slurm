#! /bin/bash

#SBATCH -J hierachy_separate_file_seqkit
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 2
#SBATCH -c 30

source activate seqkit

trim_fa=$1

trim_fa_dir=`dirname ${trim_fa}` 

seqkit grep -j 30 -n -r -p "no-split|Jump" ${trim_fa} -o ${trim_fa_dir}/`basename ${trim_fa} .fasta.gz`.moresplit.fasta.gz &

seqkit grep -v -j 30 -n -r -p "no-split|Jump" ${trim_fa} -o ${trim_fa_dir}/`basename ${trim_fa} .fasta.gz`.nomoresplit.fasta.gz &

wait
