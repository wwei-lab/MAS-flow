#! /bin/bash

#SBATCH -J combine_match
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 2
#SBATCH -c 2

source activate jupyter

python scripts/combine_alncid_match.v2.py 
