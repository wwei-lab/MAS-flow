#! /bin/bash

#SBATCH -J classify
#SBATCH -p q1
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --mem=200G
#SBATCH --ntasks 80

#Rscript scripts/classify.R
Rscript scripts/segment-classify.R
