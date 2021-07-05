#!/bin/bash
#SBATCH -N 1
#SBATCH -n 6

module add perl/5.24.1
module load vcftools/0.1.13

for i in {1..140}; do
	cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/4_merge_and_annotate_vcf_files
	u_var="$(head -n $i samp_iter | tail -n 1)"
	cd /net/gelbart/data/jacks/zacks_area/snp_mafs_and_vcfs/vcfs/
	#reorder the snv VCF columns to match the order of template.vcf
	#like this: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
	vcf-shuffle-cols -t template.vcf reform_"$u_var".vcf > temp.vcf	
	cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/4_merge_and_annotate_vcf_files/1_union_indel_vcfs
	cat "$u_var"_union.vcf > temp2.vcf
	#combine the SNV variants with the intersected indel PASS variants
	#the grep command also removes a subset of non-PASS SNV variants. Other non-PASS SNVs are filtered in a later script
	vcf-concat /net/gelbart/data/jacks/zacks_area/snp_mafs_and_vcfs/vcfs/temp.vcf /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/4_merge_and_annotate_vcf_files/1_union_indel_vcfs/temp2.vcf | grep -v -E "t_lod_fstar|alt_allele_in_normal" > /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/4_merge_and_annotate_vcf_files/2_merged_union_indel_snvs_vcfs/"$u_var"_merged.vcf
	rm /net/gelbart/data/jacks/zacks_area/snp_mafs_and_vcfs/vcfs/temp.vcf temp2.vcf ;
done


