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
genepy_df = pd.read_table(genepy_matrix_file, index_col = 'Sampleid')
    
# normalise genepy across genes using Min max scale, every value to be between 0 and 1.
def norm_genepy():
    '''
    Function to normalise down genes 
    '''
    return genepy_df.apply(preprocessing.minmax_scale, 0)
    
genepy_normalized = norm_genepy()

# binerise series (required for binerise_genepy)
def binerise_series(X, percent=95):
    '''
    Function to binerise series.
    '''
    arr1 = []
    percentile = np.percentile(X,percent)
    for row in X:
        if row > percentile:
            arr1.append(1)
        else:
            arr1.append(0)
    return arr1

# binerise scores 
def binerise_genepy(percent=95):
    '''
    Function to binerise genepy scores.
    '''
    return genepy_df.apply(binerise_series,0, percent=percent)

genepy_bin_95 = binerise_genepy()
    

# read LOEUF scores into df
LOEUF_df = pd.read_excel(loeuf_excel_file, 
                         sheet_name="Sheet1")

# create series and dict of loeuf scores upper
LOEUF_series = pd.Series(LOEUF_df["oe_lof_upper"].values, 
                         index = LOEUF_df["gene"])

LOEUF_dict = LOEUF_series.to_dict()


# removes 'nan' values
#LOEUF_dict = {k:v if not np.isnan(v) else 1 for k,v in
#              LOEUF_dict.items() }

def loeuf_score(df):
    list_LOEUF_scores = []
    # check that gene name in df is in the LOEUF dict
    for i in genepy_df.columns:
        if i in LOEUF_dict:
            # if it is, append the dictionary value
            list_LOEUF_scores.append(LOEUF_dict[i])
        else:
            # if not, append NaN
            list_LOEUF_scores.append(np.nan)
            
    # change list back into series with df.column names     
    LOEUF_score_series = pd.Series(list_LOEUF_scores, index = df.columns)
    # Fill any NaN values so that it can be devided 
    filled_series = LOEUF_score_series.fillna(1)
    return filled_series

'''
# gets a list of loeuf scores for genepy matrix
def get_loeuf_score(df_for_loeuf):
'''
    # function to get loeuf scores for df containing gene names as column headings
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
'''
    
# weight normalised scores by LOEUF
def divide_df_by_loeuf(df):
    '''
    function to divide df rows by loeuf score for each gene
    '''
    return df/loeuf_score(df)
    
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
