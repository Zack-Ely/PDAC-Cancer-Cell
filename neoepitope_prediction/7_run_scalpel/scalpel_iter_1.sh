#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8

module load python/2.7.13

#note that iteration_file_for_mutation_calling.txt contains sample IDs and file directories

for i in {1..9}; do 
	cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/1_scalpel/tcga
	var_1=$(cut -d " " -f 3 iteration_file_for_mutation_calling.txt | cut -d "-" -f 3 | sort | uniq -d | head -n $i | tail -n 1)
	#get tumor sample ID (usually ending in 01A)
        var_tum=$(grep $var_1 iteration_file_for_mutation_calling.txt | grep "\-01A")
	#get normal sample ID (usually not ending in 01A and often ending in 10A)
        var_norm=$(grep $var_1 iteration_file_for_mutation_calling.txt | grep -v "\-01A")
        tum_dir=$(echo $var_tum | cut -d " " -f 2 | sed 's/"//g')
        norm_dir=$(echo $var_norm | cut -d " " -f 2 | sed 's/"//g')
        cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/1_scalpel/tcga/results_scalpel
        mkdir $var_1
	#note the scalpel command uses a bed file derived from CGHub bitbucket - see manuscript for detail
	scalpel-discovery --somatic \
	--tumor /net/gelbart/data/jacks/TCGA/TCGA_PAAD_Exomes_2_2020/"$tum_dir"/*bam \
	--normal /net/gelbart/data/jacks/TCGA/TCGA_PAAD_Exomes_2_2020/"$norm_dir"/*bam \
	--bed /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/1_scalpel/hglft_genome_5323d_9cd6a0_final_new.bed \
	--ref /net/gelbart/data/jacks/zacks_area/gdc_ref/GRCh38.d1.vd1.fa \
	--dir ./"$var_1" \
	--two-pass ;
done
