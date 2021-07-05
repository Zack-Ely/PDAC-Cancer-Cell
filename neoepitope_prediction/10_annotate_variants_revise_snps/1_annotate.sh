#!/bin/bash
#SBATCH -N 1
#SBATCH -n 6

module add perl/5.24.1
export PATH="/net/gelbart/data/jacks/zacks_area/tools/tabix-0.2.6:$PATH";

#run this script on your VCF files containing merged SNVs and indels. This tool will annotate all of the variants.
#refer to Ensembl VEP documentation to obtain the appropriate plugins and caches. pVACTools website is also helpful

for i in {1..140} ; do
        cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/4_merge_and_annotate_vcf_files/2_merged_union_indel_snvs_vcfs
        x_var="$(head -n $i samp_iter | tail -n 1)"
	/net/gelbart/data/jacks/zacks_area/tools/ensembl-vep/vep --input_file /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/4_merge_and_annotate_vcf_files/2_merged_union_indel_snvs_vcfs/"$x_var"_merged.vcf --output_file "$x_var"_annotated.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta /net/gelbart/data/jacks/zacks_area/gdc_ref/GRCh38.d1.vd1.fa --offline --cache --plugin Downstream --plugin Wildtype --dir_cache /net/gelbart/data/jacks/zacks_area/ensemble_vep_cache --dir_plugins /net/gelbart/data/jacks/zacks_area/5_apr_att/VEP_plugins --transcript_version --pick --pick_order rank,canonical,appris,tsl,biotype,ccds,length,mane --force_overwrite ;
done

