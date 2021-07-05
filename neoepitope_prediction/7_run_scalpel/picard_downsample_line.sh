java -jar $PICARD DownsampleSam I=/net/gelbart/data/jacks/TCGA/TCGA_PAAD_Exomes_2_2020/"$tumor_var"/"$tumor_bam" O=/net/gelbart/data/jacks/zacks_area/downsample_"$patient_id".bam P=0.5
