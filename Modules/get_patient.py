import pandas as pd
import numpy as np


# get data from get_data.py
from get_data import genepy_norm_loeuf

# function to get genes and scores for a patient
def get_patient_as_df(df, patientID):
    '''
    function to return individual patient data, can use patientID or index number. Output is a df.
    '''
    if type(patientID) == str:
        gene_df  = df.loc[patientID].to_frame(name='gen_score')
    else:
        gene_df = df.iloc[patientID].to_frame(name='gen_score')
    # creates an extra column
    # x['gene'] = x.index
    return gene_df


# function to get genes and scores for a patient
def get_patient_as_Series(df, patientID):
    '''
    function to return individual patient data, can use patientID or index number. Output is a series.
    '''
    if type(patientID) == str:     # checks if patientID input is string or integer
        patient_series = df.loc[patientID]
    else:
        patient_series = df.iloc[patientID]
    # creates an extra column
    # x['gene'] = x.index
    return patient_series



# return scores for patient above 0
def get_scores_above_zero(df, patientID):
    '''
    function to return scores above 0, can use patientID or index number.
    '''
    if type(patientID) == str:
        for i in df.loc[patientID]:
            if i > 0:
                print(i)
    else:
        for i in df.iloc[patientID]:
            if i > 0:
                print(i)

## return series for patient data greater than 0
def get_scores_above_zero_as_series(df, patientID):
    '''
    function to return series above 0, can use patientID or index number. Output is a series.
    '''
    patient_series_bool = get_patient_as_Series(df, patientID) > 0
    patient_series_grt_zero = get_patient_as_Series(df, patientID).iloc[np.where(patient_series_bool)[0]]
    return patient_series_grt_zero


# get gene names greater than 0 for patient
def get_genes_above_zero(df, patientID):
    '''
    Function to return gene names with a score > 0 for a selected patient. Can use patientID or index number. Output is an index array. 
    '''
    return (get_scores_above_zero_as_series(df, patientID)
            .index
           )


# get top 2000 gene names
def get_top_2000_genes(df, patientID):
    '''
    Function to sort and return top 2000 gene names for a selected patient. 
    '''
    return (get_patient_as_Series(df, patientID)
            .sort_values(ascending=False)[:2000]
            .index
           )



# get string IDs
# TODO add no_of_proteins variable + test - link to get_PPI
def get_string_IDs(df, patientID, no_of_proteins = 2000):
    '''
    Function to return the StringDB IDs. Makes for quicker pull down of PPIs.
    Note: When using the StringDB API to pull down PPI interactions, max query of proteins is 2000, hence max gene names is 2000.  
    '''
    
    import requests

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"

    # import gene names
    gene_list = get_top_2000_genes(df, patientID)

    # human NCBI identifier 9606
    species_id = 9606

    ## Set parameters
    params = {
        "identifiers" : "\r".join(gene_list), # your protein list
        "species" : species_id, # species NCBI identifier 
        "limit" : 1, # only one (best) identifier per input protein
        "echo_query" : 1, # see your input identifiers in the output
        "caller_identity" : "Research_Project" # your app name
    }

    ## Construct URL
    request_url = "/".join([string_api_url, output_format, method])

    ## Call STRING
    results = requests.post(request_url, data=params)
    string_id = []
    for line in results.text.strip().split("\n"):
        l = line.split("\t")
        string_identifier = l[2]
        string_id.append(string_identifier)

    return string_id[:no_of_proteins]



def get_patient_phenotype(patientID):
    '''
    function that retrieves patient phenotype. Returns as df.  
    '''
    from get_data import patient_phenotype
    return patient_phenotype.loc[[patientID]]


def get_gene_phenotype(df, gene):
    '''
    function that retireves patient phenotype data from get_data.py and returns the phenotypes accosiated for that given gene.
    '''
    from get_data import patient_phenotype as patient_phe
    
    pheno_df_hold = pd.DataFrame(columns=['gender','Diagnosis'])
    gene_df = df[gene].sort_values(ascending=False).to_frame()
    
    for patient_x in gene_df.index:
        for patient_y in patient_phe.index:
            if patient_x == patient_y:
                s = patient_phe.loc[patient_x]
                pheno_df_hold = pheno_df_hold.append(s)
            else:
                continue
    return(gene_df.join(pheno_df_hold))
