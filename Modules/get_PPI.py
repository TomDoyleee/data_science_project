import pandas as pd
# import get string ids
import get_patient


def get_protein_interactions(df, patientID, no_of_proteins = 2000):
    '''
    function returns the StringDB API response for PPI data. 
    '''
    import requests
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    ## Construct URL
    request_url = "/".join([string_api_url, output_format, method])

    # set confience_score between 0-1
    species_id = 9606

    ## Set parameters
    my_genes = get_patient.get_string_IDs(df, patientID, no_of_proteins = 2000)
    
    params = {
        "identifiers" : "%0d".join(my_genes), # your protein
        "species" : species_id, # species NCBI identifier 
        "caller_identity" : "Research_Project" # your app name
    }

    """ other params
"required_score" :	# threshold of significance to include a interaction, a number between 0 and 1000 (default depends on the network)
    "network_type" :	# network type: functional (default), physical
    "add_nodes"	: # adds a number of proteins with to the network based on their confidence score
    "show_query_node_labels" :	# when available use submitted names in the preferredName column when (0 or 1) (default:0)
    """

    ## Call STRING
    response = requests.post(request_url, data=params)
    return response


def get_PPI_df(df, patientID, no_of_proteins = 2000):
    '''
    Function pulls down PPI data from StringDB. Output is a df.
    '''
    response = get_protein_interactions(df, patientID, no_of_proteins = 2000)
    patient_df = pd.DataFrame(columns = ['stringId_A',
                                         'stringId_B',
                                         'preferredName_A', 
                                         'preferredName_B',
                                         'ncbiTaxonId',
                                         'score',
                                         'nscore',
                                         'fscore',
                                         'pscore',
                                         'ascore',
                                         'escore',
                                         'dscore', 
                                         'tscore'])
  
    
    for line in response.text.strip().split("\n"):
        # seperates each line into list
        l = line.strip().split("\t")

        # adds each line to the data frame 
        patient_df.loc[len(patient_df)] = l
        
    patient_df = patient_df.drop_duplicates(subset=['preferredName_A',
                                                    'preferredName_B'],
                                            keep='first')
    return patient_df



def get_PPI_with_confidence(df, patientID, confidence_score = 0.6):
    '''
    Function to return PPI data to specific confidence score. 
    '''
    response = get_protein_interactions(df, patientID)
    patient_df = pd.DataFrame(columns = ['stringId_A', 
                                         'stringId_B',
                                         'preferredName_A',
                                         'preferredName_B',
                                         'ncbiTaxonId',
                                         'score',
                                         'nscore',
                                         'fscore',
                                         'pscore',
                                         'ascore',
                                         'escore',
                                         'dscore',
                                         'tscore'])
    
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        
        ## filter the interaction according to combined score
        combined_score = float(l[5])
        
        if combined_score > confidence_score:
            patient_df.loc[len(patient_df)] = l
    patient_df = patient_df.drop_duplicates(subset=['preferredName_A',
                                                    'preferredName_B'],
                                            keep='first')
    return(patient_df)

