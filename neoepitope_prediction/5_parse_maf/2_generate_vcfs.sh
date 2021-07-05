#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8

module add perl/5.24.1

#this script will convert patient maf files to vcf file format. It uses a perl script available in a github repository (see manuscript)
#refer to the tool's site for more information

for i in {1..36} ; do
	cd /net/gelbart/data/jacks/zacks_area/second_pdac_cohort/5_parse_mutect_maf_extended/TCGA/snp_mafs_and_vcfs
	m_var="$(ls *maf | head -n $i | tail -n 1)"
	perl /net/gelbart/data/jacks/zacks_area/tools/mskcc-vcf2maf-5453f80/maf2vcf.pl --input-maf "$m_var" --output-dir  /net/gelbart/data/jacks/zacks_area/second_pdac_cohort/5_parse_mutect_maf_extended/TCGA/snp_mafs_and_vcfs/vcfs --per-tn-vcfs --ref-fasta /net/gelbart/data/jacks/zacks_area/gdc_ref/GRCh38.d1.vd1.fa ;
done
