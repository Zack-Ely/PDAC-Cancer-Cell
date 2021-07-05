#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8

conda init bash
module load miniconda2/4.7.12
source activate pvactool

for i in {1..3} ; do
	#iterate over a list of sample IDs and set variable
	u_var="$(head -n $i /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/5_predict_neoantigens/package_uger/list_all | tail -n 1)"
	#Change the column name 'TUMOR' of VCF file to sample ID for pVACTools input
	sed "s/TUMOR/$u_var/g" /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/5_predict_neoantigens/package_uger/revised_"$u_var"_annotated.vcf > /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/5_predict_neoantigens/package_uger/"$u_var"_temp.vcf
	#note filters are included below but pVACtools still outputs a file containing all raw results, which we filter after the commands complete for all samples
	pvacseq run -e 8,9,10,11 --iedb-install-directory /net/gelbart/data/jacks/zacks_area/ -b 1000 -m median -k --net-chop-method cterm --netmhc-stab -a sample_name -d full \
	--normal-sample-name NORMAL -c 1.0 --tdna-cov 5 --tdna-vaf 0.1 /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/5_predict_neoantigens/package_uger/"$u_var"_temp.vcf "$u_var" $(grep "$u_var" /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/5_predict_neoantigens/package_uger/ultimate_master_matrix_generation_all.csv | cut -d "," -f 2-9 | sed 's/A/HLA-A/g' | sed 's/B/HLA-B/g' | sed 's/C/HLA-C/g' | sed 's/E/HLA-E/g') NetMHC NetMHCpan SMM SMMPMBEC /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/5_predict_neoantigens/trial_outs/"$u_var"/
	rm /net/gelbart/data/jacks/zacks_area/third_pdac_analysis/5_predict_neoantigens/package_uger/"$u_var"_temp.vcf ;
done
