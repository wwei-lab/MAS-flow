#! /bin/bash

#SBATCH -J minimap2_trim
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 1
#SBATCH --mem=800G
#SBATCH -c 100

source activate minimap2

fa=$1
out_dir=$2
genome=$3

minimap2 -K500m -t 100 --secondary=no -a -x splice --splice-flank=yes ${genome} $fa |samtools sort -@ 80 > ${out_dir}/`basename $fa`.sorted.bam
