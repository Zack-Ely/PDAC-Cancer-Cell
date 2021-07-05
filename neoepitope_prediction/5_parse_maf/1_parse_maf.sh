#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8

#example script to parse the mutect MAF file obtained from TCGA for all patients into patient-specific files
#note this requires a file, list_of_final_tcga, that contains a portion of the sample ID

cd /net/gelbart/data/jacks/zacks_area/second_pdac_cohort/5_parse_mutect_maf_extended/TCGA

rm heade

#take header lines from the maf file
head -n 6 TCGA.PAAD.mutect.DOWNLOAD.somatic.maf > heade

for i in {1..36}; do
	u_var="$(head -n $i list_of_final_tcga | tail -n 1)"
	grep "$u_var"-01A TCGA.PAAD.mutect.DOWNLOAD.somatic.maf > taile
	#take only SNPs from the MAF file
	grep -v "#" taile | awk '$10 == "SNP"' > taile2
	cat heade taile2 > /net/gelbart/data/jacks/zacks_area/second_pdac_cohort/5_parse_mutect_maf_extended/TCGA/snp_mafs_and_vcfs/"$u_var".maf
	rm taile taile2 ;
done

