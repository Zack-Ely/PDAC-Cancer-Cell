#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8

module load vcftools/0.1.13

cd /net/gelbart/data/jacks/zacks_area/

#this script parses the raw strelka results by subsetting PASSed variants and outputting them to a new directory for all patients

for i in {1..28}; do
	cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/2_strelka/
	u_var="$(head -n $i iter_parse_tcga.txt | tail -n 1)"
	cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/2_strelka/tcga/results_strelka/"$u_var"/results/variants/
	vcftools --gzvcf somatic.indels.vcf.gz --remove-filtered-all --recode --recode-INFO-all --out /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/3_strelka_scalpel_mutation_processing/indel_vcfs/out_"$u_var" --stdout | gzip -c > /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/3_strelka_scalpel_mutation_processing/indel_vcfs/"$u_var"_strelka.vcf.gz ;
done
