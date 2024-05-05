#! /bin/bash

#SBATCH -J count_accuracy
#SBATCH -p q1,q2,q3
#SBATCH -N 1
#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err
#SBATCH --ntasks 100

source activate jupyter

#Rscript scripts/count_cigar.R BGI_benchmark/Nanopore_T44_0116/FAX23635_pass_12b693ce_c74c2dc4.splitted.fasta.gz.sorted.bam BGI_benchmark/Nanopore_T44_0116/ 2

Rscript scripts/count_cigar.R BGI_benchmark/BGI_data/TB20003420-202401221650300_read.splitted.fasta.gz.sorted.bam BGI_benchmark/BGI_data/ 2

