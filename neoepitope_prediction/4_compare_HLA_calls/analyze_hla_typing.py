import pandas as pd

key_seq_opt_mat = pd.read_csv("/net/gelbart/data/jacks/zacks_area/second_pdac_cohort/4_hla_analysis/TCGA/FINAL_Seq2_Optitype_comparison_matrix.csv") #make input

def remove_ambiguity(input_str):
    """This function removes ambiguity marks output by seq2HLA."""
    remove_mark = "'`"
    result_str = ""
    for i in input_str:
        if i not in remove_mark:
            result_str = result_str + i
    return result_str

def eval_seq_opti_consistency(input_matr):
    """This function evaluates the consistency between seq2HLA and OptiType calls."""
    A_mismatches = 0
    A_matches = 0
    B_mismatches = 0
    B_matches = 0
    C_mismatches = 0
    C_matches = 0
    annotated_matr = input_matr
    for i in range(len(input_matr)):
        #A allele
        #heterozygous seq2hla
        if remove_ambiguity(input_matr.iloc[i, 1]) != remove_ambiguity(input_matr.iloc[i, 3]):
            if remove_ambiguity(input_matr.iloc[i, 1]) == remove_ambiguity(input_matr.iloc[i, 13]) or remove_ambiguity(input_matr.iloc[i, 1]) == remove_ambiguity(input_matr.iloc[i, 14]):
                A_matches = A_matches + 1
            else:
                A_mismatches = A_mismatches + 1
                annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon' #add a label to mark inconsistent calls
            if remove_ambiguity(input_matr.iloc[i, 3]) == remove_ambiguity(input_matr.iloc[i, 13]) or remove_ambiguity(input_matr.iloc[i, 3]) == remove_ambiguity(input_matr.iloc[i, 14]):
                A_matches = A_matches + 1
            else:
                A_mismatches = A_mismatches + 1
                annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
        #homozygous seq2hla
        else:
            #if optitype is also homozygous
            if remove_ambiguity(input_matr.iloc[i, 13]) == remove_ambiguity(input_matr.iloc[i, 14]):
                #add 2 matches if all alleles match
                if remove_ambiguity(input_matr.iloc[i, 1]) == remove_ambiguity(input_matr.iloc[i, 13]) or remove_ambiguity(input_matr.iloc[i, 1]) == remove_ambiguity(input_matr.iloc[i, 14]):
                    A_matches = A_matches + 2
                else:
                    A_mismatches = A_mismatches + 2
                    annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
            #optitype is heterozygous
            else:
                #both heterozygous opti alleles do not match seq2 homozygous allele
                if (not(remove_ambiguity(input_matr.iloc[i, 13]) == (remove_ambiguity(input_matr.iloc[i, 1]))) and not(remove_ambiguity(input_matr.iloc[i, 14]) == (remove_ambiguity(input_matr.iloc[i, 3])))):
                    A_mismatches = A_mismatches + 2
                    annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
                else:
                    A_mismatches = A_mismatches + 1
                    A_matches = A_matches + 1
                    annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
        #B allele
        #heterozygous seq2hla
        if remove_ambiguity(input_matr.iloc[i, 5]) != remove_ambiguity(input_matr.iloc[i, 7]):
            if remove_ambiguity(input_matr.iloc[i, 5]) == remove_ambiguity(input_matr.iloc[i, 15]) or remove_ambiguity(input_matr.iloc[i, 5]) == remove_ambiguity(input_matr.iloc[i, 16]):
                B_matches = B_matches + 1
            else:
                B_mismatches = B_mismatches + 1
                annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
            if remove_ambiguity(input_matr.iloc[i, 7]) == remove_ambiguity(input_matr.iloc[i, 15]) or remove_ambiguity(input_matr.iloc[i, 7]) == remove_ambiguity(input_matr.iloc[i, 16]):
                B_matches = B_matches + 1
            else:
                B_mismatches = B_mismatches + 1
                annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
        #homozygous seq2hla
        else:
            #if optitype is also homozygous
            if remove_ambiguity(input_matr.iloc[i, 15]) == remove_ambiguity(input_matr.iloc[i, 16]):
                #add 2 matches if all alleles match
                if (remove_ambiguity(input_matr.iloc[i, 5]) == remove_ambiguity(input_matr.iloc[i, 15]) or remove_ambiguity(input_matr.iloc[i, 5]) == remove_ambiguity(input_matr.iloc[i, 16])):
                    B_matches = B_matches + 2
                else:
                    B_mismatches = B_mismatches + 2
                    annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
            #optitype is heterozygous
            else:
                #both heterozygous opti alleles do not match seq2 homozygous allele
                if (not(remove_ambiguity(input_matr.iloc[i, 15]) == (remove_ambiguity(input_matr.iloc[i, 5]))) and not(remove_ambiguity(input_matr.iloc[i, 16]) == (remove_ambiguity(input_matr.iloc[i, 5])))):
                    B_mismatches = B_mismatches + 2
                    annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
                else:
                    B_mismatches = B_mismatches + 1
                    B_matches = B_matches + 1
                    annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
        #C allele
        #heterozygous seq2hla
        if remove_ambiguity(input_matr.iloc[i, 9]) != remove_ambiguity(input_matr.iloc[i, 11]):
            if remove_ambiguity(input_matr.iloc[i, 9]) == remove_ambiguity(input_matr.iloc[i, 17]) or remove_ambiguity(input_matr.iloc[i, 9]) == remove_ambiguity(input_matr.iloc[i, 18]):
                C_matches = C_matches + 1
            else:
                C_mismatches = C_mismatches + 1
                annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
            if remove_ambiguity(input_matr.iloc[i, 11]) == remove_ambiguity(input_matr.iloc[i, 17]) or remove_ambiguity(input_matr.iloc[i, 11]) == remove_ambiguity(input_matr.iloc[i, 18]):
                C_matches = C_matches + 1
            else:
                C_mismatches = C_mismatches + 1
                annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
        #homozygous seq2hla
        else:
            #if optiypte is also homozygous
            if remove_ambiguity(input_matr.iloc[i, 17]) == remove_ambiguity(input_matr.iloc[i, 18]):
                #add 2 matches if all alleles match
                if (remove_ambiguity(input_matr.iloc[i, 9]) == remove_ambiguity(input_matr.iloc[i, 17]) or remove_ambiguity(input_matr.iloc[i, 9]) == remove_ambiguity(input_matr.iloc[i, 18])):
                    C_matches = C_matches + 2
                else:
                    C_mismatches = C_mismatches + 2
                    annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
            #optitype is heterozygous
            else:
                #both heterozygous opti alleles do not match seq2 homozygous allele
                if (not(remove_ambiguity(input_matr.iloc[i, 17]) == (remove_ambiguity(input_matr.iloc[i, 9]))) and not(remove_ambiguity(input_matr.iloc[i, 18]) == (remove_ambiguity(input_matr.iloc[i, 9])))):
                    C_mismatches = C_mismatches + 2
                    annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
                else:
                    C_mismatches = C_mismatches + 1
                    C_matches = C_matches + 1
                    annotated_matr.iloc[i, 0] = annotated_matr.iloc[i, 0] + '_incon'
    return A_matches, A_mismatches, B_matches, B_mismatches, C_matches, C_mismatches

def resolve_inconsistencies(inp_matr):
    """replace seq2HLA calls with Optitype calls where there are inconsistencies."""
    working_matrix = inp_matr
    for i in range(len(working_matrix)):
        if 'incon' in working_matrix.iloc[i, 0]:
            working_matrix.iloc[i, 1] = working_matrix.iloc[i, 13]
            working_matrix.iloc[i, 3] = working_matrix.iloc[i, 14]
            working_matrix.iloc[i, 5] = working_matrix.iloc[i, 15]
            working_matrix.iloc[i, 7] = working_matrix.iloc[i, 16]
            working_matrix.iloc[i, 9] = working_matrix.iloc[i, 17]
            working_matrix.iloc[i, 11] = working_matrix.iloc[i, 18]
    working_matrix = working_matrix.iloc[:, [0,1,3,5,7,9,11]]
    return working_matrix

#below quantifies mismatches and implicitly adds 'incon' to the key_seq_opt_matrix
consistency = eval_seq_opti_consistency(key_seq_opt_mat)
filter_matrix = resolve_inconsistencies(key_seq_opt_mat)
filter_matrix.to_csv('/net/gelbart/data/jacks/zacks_area/second_pdac_cohort/4_hla_analysis/TCGA/final_filtered_hla_matrix.csv', index=False)
print(consistency)
#now you can examine the result matrix for the number of instances of "incon" to determine consistency
