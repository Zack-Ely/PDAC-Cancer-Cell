#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16

#Note: to run on the PancSeq cohort, make parameter modifications as described in our manuscript

module add r/3.5.1
module load python/2.7.13

for i in {1..10}; do
        cd /net/gelbart/data/jacks/TCGA/TCGA_PAAD_Exomes_2_2020
        cp ./$(ls | head -n $i | tail -n 1)/*bam /net/gelbart/data/jacks/zacks_area/
        cp ./$(ls | head -n $i | tail -n 1)/*bai /net/gelbart/data/jacks/zacks_area/
	cd ./$(ls | head -n $i | tail -n 1)/
        l_var="$(ls *bam)"
	j_var="$(ls *bai)"
	g_var="$(echo $l_var | cut -d "." -f 2)"
        cd /net/gelbart/data/jacks/zacks_area/
	#prepare bam files for input to OptiType (see our manuscript and OptiType documention)
	samtools view -b -X "$l_var" "$j_var" chr6:28510120-33480577 -o 2_"$l_var"
	samtools sort -n 2_"$l_var" -o 2_"$l_var".qsort
	#convert sorted bam file to fastq files
	bedtools bamtofastq -i 2_"$l_var".qsort -fq ./2_"$l_var".end1.fq -fq2 ./2_"$l_var".end2.fq
	razers3 -i 95 -m 1 -dr 0 -o "$l_var"_fished_1.bam /net/gelbart/data/jacks/zacks_area/tools/OptiType/data/hla_reference_dna.fasta 2_"$l_var".end1.fq
	samtools bam2fq "$l_var"_fished_1.bam > 1_"$l_var"_fished.fastq
	razers3 -i 95 -m 1 -dr 0 -o "$l_var"_fished_2.bam /net/gelbart/data/jacks/zacks_area/tools/OptiType/data/hla_reference_dna.fasta 2_"$l_var".end2.fq
	samtools bam2fq "$l_var"_fished_2.bam > 2_"$l_var"_fished.fastq
	mkdir /net/gelbart/data/jacks/zacks_area/results_optitype/out_"$g_var"
	#note the python script used below is available from OptiType
	python /net/gelbart/data/jacks/zacks_area/tools/OptiType/OptiTypePipeline.py -i 1_"$l_var"_fished.fastq 2_"$l_var"_fished.fastq --dna -c /net/gelbart/data/jacks/zacks_area/tools/OptiType/config.ini --outdir /net/gelbart/data/jacks/zacks_area/results_optitype/out_"$g_var"
	rm 2_"$l_var"_fished.fastq 1_"$l_var"_fished.fastq "$l_var"_fished_1.bam "$l_var"_fished_2.bam 2_"$l_var" 2_"$l_var".end2.fq 2_"$l_var".end1.fq "$l_var" "$j_var" 2_"$l_var".qsort ;
done
