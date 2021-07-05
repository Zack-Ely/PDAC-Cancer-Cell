#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8


#Manta is supposed to be run and then piped into the strelka workflow

module load python/2.7.13

for i in {1..9}; do
	cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/2_strelka/tcga
	var_1=$(cut -d " " -f 3 iteration_file_for_mutation_calling.txt | cut -d "-" -f 3 | sort | uniq -d | head -n $i | tail -n 1)
        var_tum=$(grep $var_1 iteration_file_for_mutation_calling.txt | grep "\-01A")
        var_norm=$(grep $var_1 iteration_file_for_mutation_calling.txt | grep -v "\-01A")
        tum_dir=$(echo $var_tum | cut -d " " -f 2 | sed 's/"//g')
        norm_dir=$(echo $var_norm | cut -d " " -f 2 | sed 's/"//g')
        cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/2_strelka/tcga/results_strelka
        mkdir $var_1
	#note this step requires manta results as an input parameter to Strelka.
	python /net/gelbart/data/jacks/zacks_area/tools/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
	--normalBam /net/gelbart/data/jacks/TCGA/TCGA_PAAD_Exomes_2_2020/$norm_dir/*bam \
	--tumorBam /net/gelbart/data/jacks/TCGA/TCGA_PAAD_Exomes_2_2020/$tum_dir/*bam \
	--referenceFasta /net/gelbart/data/jacks/zacks_area/gdc_ref/GRCh38.d1.vd1.fa \
	--indelCandidates /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/0_manta/tcga/results_manta/$var_1/results/variants/candidateSmallIndels.vcf.gz \
	--runDir /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/2_strelka/tcga/results_strelka/$var_1/ \
	--exome
	python /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/2_strelka/tcga/results_strelka/$var_1/runWorkflow.py -m local -j 8 ;
done
