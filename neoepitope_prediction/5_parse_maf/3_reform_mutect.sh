#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=END
#SBATCH --mail-user=zack_ely@mit.edu

#to use this script, add the list of patients to the /vcf/ directory (eg list_of_final_tcga)

#this script reformats the converted VCF files to contain a different column-header line

#note this script uses a file called gen_header that contains the standard VCF header line: 
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR	NORMAL

for i in {1..36}; do
	cd /net/gelbart/data/jacks/zacks_area/second_pdac_cohort/5_parse_mutect_maf_extended/TCGA/snp_mafs_and_vcfs/vcfs
	r_var="$(head -n $i list_of_final_tcga | tail -n 1)"
	#extract most header lines from a sample's vcf file
	grep "##" "$r_var".vcf > piece_header_lines
	grep -v "#" "$r_var".vcf > content_lines
	cat piece_header_lines gen_header content_lines > reform_"$r_var".vcf
	rm piece_header_lines content_lines ;
done
