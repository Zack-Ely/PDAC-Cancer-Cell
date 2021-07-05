#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8

module load python/2.7.13

#note that the iteration_file simply contains 3 columns with the sample_ID in the last column and specific path names in
#the other columns. The most important parts are the manta calls listed within this loop.

for i in {1..9}; do
	cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/0_manta/tcga
	var_1=$(cut -d " " -f 3 iteration_file_for_mutation_calling.txt | cut -d "-" -f 3 | sort | uniq -d | head -n $i | tail -n 1)
	var_tum=$(grep $var_1 iteration_file_for_mutation_calling.txt | grep "\-01A")
	var_norm=$(grep $var_1 iteration_file_for_mutation_calling.txt | grep -v "\-01A")
	tum_dir=$(echo $var_tum | cut -d " " -f 2 | sed 's/"//g')
	norm_dir=$(echo $var_norm | cut -d " " -f 2 | sed 's/"//g')
	cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/0_manta/tcga/results_manta/
	mkdir $var_1
	python /net/gelbart/data/jacks/zacks_area/tools/manta/bin/configManta.py \
	--normalBam /net/gelbart/data/jacks/TCGA/TCGA_PAAD_Exomes_2_2020/$norm_dir/*bam \
	--tumorBam /net/gelbart/data/jacks/TCGA/TCGA_PAAD_Exomes_2_2020/$tum_dir/*bam \
	--referenceFasta /net/gelbart/data/jacks/zacks_area/gdc_ref/GRCh38.d1.vd1.fa \
	--runDir /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/0_manta/tcga/results_manta/$var_1/ \
	--exome
	python /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/0_manta/tcga/results_manta/$var_1/runWorkflow.py -j 8 ;
done
