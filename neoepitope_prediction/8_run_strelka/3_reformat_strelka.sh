#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8


#This script can help reformat strelka output files by reformatting the Strelka-specific depth statistics into
#a more common format that can be parsed by pVACtools, etc. This script assumes you're in a directory
#that only contains strelka-outputted VCF files for your cohort of interest.

for i in {1..85}; do
	r_var="$(ls *vcf | head -n $i | tail -n 1)"
	grep "#" "$r_var" > piece_header_lines
	grep -v "#" "$r_var" | cut -f 9 | sed 's/DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50/GT:AD:DP/g' > piece_col_9
	grep -v "#" "$r_var" | cut -f 1-8 > piece_col_1_8
	grep -v "#" "$r_var" | cut -f 10 > col_10
	grep -v "#" "$r_var" | cut -f 11 > col_11
	while read line ; do
		tar="$(echo $line | cut -d ":" -f 3 | cut -d "," -f 1)"
		tir="$(echo $line | cut -d ":" -f 4 | cut -d "," -f 1)"
		dp="$(($tar+$tir))"
		echo 0/1:"$tar","$tir":"$dp" >> new_col_10 ;
	done < col_10
	while read line ; do
		tar2="$(echo $line | cut -d ":" -f 3 | cut -d "," -f 1)"
		tir2="$(echo $line | cut -d ":" -f 4 | cut -d "," -f 1)"
		dp2="$(($tar2+$tir2))"
		echo 0/1:"$tar2","$tir2":"$dp2" >> new_col_11 ;
	done < col_11
	paste piece_col_1_8 piece_col_9 new_col_10 new_col_11 > reform_bottom
	cat piece_header_lines reform_bottom > reformed_"$r_var"
	rm piece_col_9 piece_col_1_8 col_10 col_11 new_col_10 new_col_11 reform_bottom piece_header_lines ;
done
