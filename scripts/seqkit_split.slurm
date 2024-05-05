#! /bin/bash

#SBATCH -J seqkit_split
#SBATCH -p q1
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mem=50G
#SBATCH --ntasks 1
#SBATCH -c 100

source activate seqkit

seqkit split -j 100 /cluster/home/liuhengxin/rawdata/ONT-15LIGA-240312/BMK240116-BV827-ZX01-0401/BMK_DATA_20240312105243_1/ONT_BV827-001N0001_clean.fq.gz -s 10000000
