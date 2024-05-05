#! /bin/bash

#cat /cluster/home/liuhengxin/rawdata/ONT-15LIGA-240312/BMK240116-BV827-ZX01-0401/BMK_DATA_20240312105243_1/ONT_BV827-001N0001_clean.fq.gz.split/fqlist.txt|while read id
cat /cluster/home/liuhengxin/rawdata/ONT-15LIGA-240312/BMK240116-BV827-ZX01-0401/BMK_DATA_20240312105243_1/ONT_BV827-001N0001_clean.fq.gz.split/fqlist.txt|while read id
do
        echo ${id}
	fastq=${id}
        fastq_base=`basename ${fastq}`
        fastq_dir=`dirname ${fastq}`
	outdir=data/`basename $id`
        mkdir -p ${outdir}
	genome=../genome/GCF_000364345.1_Macaca_fascicularis_5.0_genomic.fna
	mismatch=2
	moresplit_mismatch=3

#	trim_fa=${fastq_dir}/"${fastq_base%%.*}".trimmed.fasta.gz
	trim_fa=${fastq_dir}/`basename $fastq_base .fq.gz`.trimmed.fasta.gz
	linkerids=../data/adapterids.linker.txt
	linkerseqs=../data/adapterseqs.linker.txt

	CIDadapter_ids=../data/adapterids.head.txt
	CIDadapter_seqs=../data/adapterseqs.head.txt

	moresplit_trim_fa=`dirname ${trim_fa}`/`basename ${trim_fa} .trimmed.fasta.gz`.trimmed.moresplit.fasta.gz

	moresplit_trim_fa2=`dirname ${trim_fa}`/`basename ${moresplit_trim_fa} .trimmed.moresplit.fasta.gz`.trimmed.moresplit.moresplit.fasta.gz

	nomoresplit_trim_fa=`dirname ${trim_fa}`/`basename ${trim_fa} .trimmed.fasta.gz`.trimmed.nomoresplit.fasta.gz

	final=`dirname ${trim_fa}`/`basename ${trim_fa} .trimmed.fasta.gz`.splitted.fasta.gz

        ncpus=100

	# 01 find adapter
#	jobid1=$(sbatch scripts/find_adapter.slurm ${fastq} ${mismatch} ${outdir} ${linkerids} ${linkerseqs} | awk '{print $4}')

	# 02 segment classify
#       -d afterok:$jobid1
#	jobid2=$(sbatch scripts/segment-classify.slurm ${outdir} ${fastq} ${mismatch}| awk '{print $4}')

	# 03 do splitting
	# -d afterok:$jobid2
	#jobid3=$(sbatch  scripts/do_splitting.slurm ${outdir} ${fastq} ${mismatch}| awk '{print $4}')

	# 04 moresplit seqkit
	#jobid4=$(sbatch -d afterok:$jobid3  scripts/hierachy_separate_file.slurm ${trim_fa}| awk '{print $4}')

	# 05 moresplit find adapter
	#jobid5=$(sbatch -d afterok:$jobid4  scripts/find_adapter.v2.slurm ${moresplit_trim_fa} ${moresplit_mismatch} ${outdir} ${linkerids} ${linkerseqs}| awk '{print $4}')

	# 06 moresplit segment classify
	#jobid6=$(sbatch -d afterok:$jobid5  scripts/segment-classify.v2.slurm ${moresplit_mismatch} ${outdir}| awk '{print $4}')

	# 07 moresplit segment splitting
	#jobid7=$(sbatch -d afterok:$jobid6 scripts/do_splitting.v2.slurm ${outdir} ${moresplit_trim_fa} ${moresplit_mismatch}| awk '{print $4}')

	# 08 merge
	#jobid8=$(sbatch -d afterok:$jobid7 scripts/merge_fasta.slurm ${moresplit_trim_fa2} ${nomoresplit_trim_fa} ${trim_fa}| awk '{print $4}')

	# 09 minimap2
	#joibid9=$(sbatch -d afterok:$jobid8 scripts/minimap2.slurm ${final} ${outdir} ${genome})

	# 10_1 extract CID-find CID oligo
        # -d afterok:$jobid9
	#jobid10=$(sbatch scripts/find_adapter.v3.slurm ${final} ${mismatch} ${outdir} ${CIDadapter_ids} ${CIDadapter_seqs}  | awk '{print $4}')

	# 10_2 extract CID-extract CID sequence
        #-d afterok:$jobid10
	$(sbatch scripts/extractCID.slurm ${outdir} ${mismatch} ${ncpus} ${final} ${CIDadapter_ids}  CIDcapture CIDoligo)

done
