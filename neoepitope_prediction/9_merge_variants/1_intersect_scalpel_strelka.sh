#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8

cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/4_merge_and_annotate_vcf_files

while read line; do 

	cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/3_strelka_scalpel_mutation_processing/indel_vcfs
	#take the 1st 5 columns of the variant lines from strelka file
	grep -v "#" reformed_"$line"_strelka.vcf | cut -f 1-5 > temp
	#take all of the variant lines from strelka file
	grep -v "#" reformed_"$line"_strelka.vcf > temp2
	#take all header lines
	grep "#" reformed_"$line"_strelka.vcf > temp3
	#moves these temporary files into the scalpel results directory
	mv ./temp /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/1_scalpel/tcga/results_scalpel/"$line"/

	mv ./temp2 /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/1_scalpel/tcga/results_scalpel/"$line"/

	mv ./temp3 /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/1_scalpel/tcga/results_scalpel/"$line"/

	cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/1_scalpel/tcga/results_scalpel/"$line"
	#use the first 5 columns from passed strelka variants to idenify variants also called PASS by Scalpel
	grep -f temp somatic.indel.vcf | grep "PASS" | cut -f 1-5 > final1
	#take the latter 6 columns of these intersected PASS variants from the original strelka file
	grep -f final1 temp2 | cut -f 6-11 > second_column_set
	#paste the first 5 and the last 6 columns for all of the mutual PASS variants
	paste final1 second_column_set > final2
	#add the strelka header back to these variants, now forming the intersected variant VCF file
	cat temp3 final2 > /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/4_merge_and_annotate_vcf_files/union_indel_vcfs/"$line"_union.vcf

	rm final1 temp temp2 temp3 final2 second_column_set;

done < samp_iter
#samp_iter is just a line-by-line list of sample IDs
