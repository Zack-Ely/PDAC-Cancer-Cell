#!/bin/bash

module add perl/5.24.1
module load vcftools/0.1.13

#Refer to methods section in manuscript for context. This script is meant to subset on PASSed SNVs, but it also
#will retain certain non-PASS designations that occur in genes associated with pancreatic cancer but which
#variants are more likely truly somatic. For example, many KRAS G12 variants were not marked PASS even though
#these are almost certainly true somatic variants driving the corresponding tumor sample. So, we retain
#such instances to help reduce recurrent false negatives

for i in {1..140} ; do
        cd /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/4_merge_and_annotate_vcf_files/2_merged_union_indel_snvs_vcfs
        u_var="$(head -n $i samp_iter | tail -n 1)"
	cat "$u_var"_annotated.vcf > temp.vcf
	grep "#" temp.vcf > top
	grep -v "#" temp.vcf | grep "PASS" | grep -v -E "KRAS|TP53|GNAS|RNF43|PLEC|FLG|AHNAK|APOB|CSMD1|PLXNA1|MCM6|MKI67|SIPA1" > m1
	grep "KRAS|" temp.vcf | grep -E "PASS|clustered_event|homologous_mapping_event|panel_of_norm" > m2
	grep TP53 temp.vcf | grep -E "PASS|panel_of_norm|bSeq" > m3
	grep "GNAS|" temp.vcf | grep -E "PASS|clustered_events|panel_of_norm" > m4
	grep RNF43 temp.vcf | grep -E "PASS|clustered_events" > m5
	grep "PLEC|" temp.vcf | grep -E "PASS|clustered_events|panel_of_norm" > m6
	grep "FLG|" temp.vcf | grep -E "PASS|panel_of_norm" > m7
	grep "AHNAK|" temp.vcf | grep -E "PASS|panel_of_norm" > m8
	grep "APOB|" temp.vcf | grep -E "PASS|panel_of_norm" > m9
	grep CSMD1 temp.vcf | grep -E "PASS|panel_of_norm" > m10
	grep PLXNA1 temp.vcf | grep -E "PASS|clustered_events" > m11
	grep MCM6 temp.vcf | grep -E "PASS|clustered_events" > m12
	grep MKI67 temp.vcf | grep -E "PASS|clustered_events" > m13
	grep SIPA1 temp.vcf | grep -E "PASS|clustered_events|homologous_mapping" > m14
	rm temp.vcf
	cat top m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 > temp.vcf
	vcf-sort temp.vcf > revised_"$u_var"_annotated.vcf # this VCF file is the one used as input to pVACTools neoantigen prediction
	rm temp.vcf top m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14;
done
