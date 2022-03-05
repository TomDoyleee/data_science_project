import pandas as pd
import numpy as np
from sklearn import preprocessing
import os

# change wd
os.chdir("/Users/your_name/your/working/directory")

# genepy matrix file in relation to wd
genepy_matrix_file = "../Data/genepy.matrix"

# loeuf excel file in relation to wd
loeuf_excel_file = "../Data/gnomad_loeuf.xlsx"

# patient phenotype table
patient_phenotype_file = '../Data/tab_delim_phenotype.txt'


# import genepy matix as df. 'Sampleid' is column name for the sample ids set to be row index.
genepy_df = pd.read_table(genepy_matrix_file,
                          index_col = 'Sampleid')
    
# normalise genepy across genes using Min max scale, every value to be between 0 and 1.
def norm_genepy():
    '''
    Function to normalise across genes for each patient
    '''
    genepy_normalized_hold = pd.DataFrame()
    
    for i in range(len(genepy_df.index)):
        row = preprocessing.minmax_scale(genepy_df.iloc[i,])
        row_series = pd.Series(row, 
                               index = genepy_df.columns)
        # adds each line to the data frame 
        genepy_normalized_hold = genepy_normalized_hold.append(row_series, 
                                                               ignore_index=True)
    #adds index values (patientID) back
    genepy_normalized_hold.index = genepy_df.index
    
    return genepy_normalized_hold
    
genepy_normalized = norm_genepy()




# read LOEUF scores into df
LOEUF_df = pd.read_excel(loeuf_excel_file, 
                         sheet_name="Sheet1")

# create series and dict of loeuf scores upper
LOEUF_series = pd.Series(LOEUF_df["oe_lof_upper"].values, 
                         index = LOEUF_df["gene"])

LOEUF_dict = LOEUF_series.to_dict()


# removes 'nan' values
LOEUF_dict = {k:v if not np.isnan(v) else 1 for k,v in
              LOEUF_dict.items() }

# gets a list of loeuf scores for genepy matrix
def get_loeuf_score(df_for_loeuf):
    '''
    function to get loeuf scores for df containing gene names as column headings
    '''
    loeuf_score_list = []
    for i in df_for_loeuf.columns:
        if i in LOEUF_dict.keys():
            if LOEUF_dict[i] == str('NaN'):
                loeuf_score_list.append(1)
            else:
                loeuf_score_list.append(LOEUF_dict[i])
        else:
            loeuf_score_list.append(1)
    return(loeuf_score_list)

# weight normalised scores by LOEUF
def divide_df_by_loeuf(df):
    '''
    function to divide df rows by loeuf score for each gene
    '''
    return df/get_loeuf_score(df)
    
genepy_norm_loeuf = divide_df_by_loeuf(genepy_normalized)
    
print(genepy_norm_loeuf.head())


# read patient phenotype into df
patient_phenotype = pd.read_table(patient_phenotype_file, index_col=0)
    

# get patient diagnosis subsets
def get_diagnosis_df(df, diagnosis):
    '''
    funtion to return subset of df based on clincal diagnosis. Diagnosis = string object e.g. 'CD'. 
    '''
    def get_diagnosis(diagnosis):
        return patient_phenotype.loc[:,'Diagnosis'] == diagnosis
    index_value = patient_phenotype.iloc[np.where(get_diagnosis(diagnosis))[0]].index
    return df.loc[index_value]

# create subset df
CD_subset = get_diagnosis_df(genepy_norm_loeuf, 'CD')
UC_subset = get_diagnosis_df(genepy_norm_loeuf, 'UC')
IBDU_subset = get_diagnosis_df(genepy_norm_loeuf, 'IBDU')
NOT_IBD_subset = get_diagnosis_df(genepy_norm_loeuf, 'NOT_IBD')

