#command line to aggregate raw neoepitope predictions across all patients in the parent output directoy
#filteing for neoepitopes with a median predicted affinity < 1000 nM (column 33 in our outpus)
cat ./*/MHC_Class_I/*.all_epitopes.tsv | awk -F "\t" '$33<1000' | wc -l
