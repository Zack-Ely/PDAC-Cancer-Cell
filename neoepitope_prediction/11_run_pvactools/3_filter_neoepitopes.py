import pandas as pd

#note this script uses an aggregate csv file constructed from the raw neoeptiope results for all of the samples. 

ant_matrix = pd.read_csv("/net/gelbart/data/jacks/zacks_area/third_pdac_analysis/5_predict_neoantigens/fixed_median_1000nm_unfiltered_matrix.csv")

def filter_antigen_matrix(antigen_matrix):
    """filter the neoantigen candidates for a variety of parameters except binding affinity but including fold change affinity"""
    filtered = antigen_matrix[(antigen_matrix["cysteine_count"] <= 1)]
    filtered_2 = filtered[(filtered["Tumor DNA Depth"] >= 5) & (filtered["Tumor DNA VAF"] >= 0.07) & ((filtered["Median Fold Change"] > 1) | (filtered["Median Fold Change"].isnull()))]
    return filtered_2

filt_median_1000 = filter_antigen_matrix(ant_matrix)

filt_median_1000.to_csv('/net/gelbart/data/jacks/zacks_area/third_pdac_analysis/5_predict_neoantigens/filtered_median_1000nm_matrix.csv', index=False)
