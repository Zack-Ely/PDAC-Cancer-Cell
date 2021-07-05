#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8

module add r/3.5.1
module load python/2.7.13

for i in {1..10}; do
	cd /net/gelbart/data/jacks/WFP_downloads/DFCI_6_2020/RNAseq_6_30_2020/RNAseq/
	#iteratively obtain paths for samples' BAM files
	i_var="$(head -n $i /net/gelbart/data/jacks/zacks_area/second_pdac_cohort/list_of_dfci_pats_for_seq2 | tail -n 1)"
	l_var="$(grep "$i_var" /net/gelbart/data/jacks/WFP_downloads/DFCI_6_2020/transfer_Pancseq_paper_cohort.csv | head -n 1 | cut -d "," -f 4 | cut -d "/" -f 2)"
	cp ./"$l_var" /net/gelbart/data/jacks/zacks_area/bam_store
	cd /net/gelbart/data/jacks/zacks_area/bam_store
	#sort bam files prior to fastq conversion
	samtools sort -n $l_var -o $l_var.qsort
	#convert bam to fastq files
	bedtools bamtofastq -i $l_var.qsort -fq $l_var.end1.fq -fq2 $l_var.end2.fq 
	gzip $l_var.end1.fq
	gzip $l_var.end2.fq
	#run seq2HLA
	python /net/gelbart/data/bcc/arjun/AJ_open/HLA-sandbox/seq2HLA/seq2HLA-master/seq2HLA.py -1 $l_var.end1.fq.gz -2 $l_var.end2.fq.gz -r "/net/gelbart/data/jacks/zacks_area/second_pdac_cohort/1_seq2hla_analysis/seq2hla_results/out_$i/out_$i"
	rm ./$l_var* ;
done
