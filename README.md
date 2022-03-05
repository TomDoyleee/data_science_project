# data_science_project

Repository for data science research project, involving downstream analysis of [Genepy](https://github.com/UoS-HGIG/GenePy-1.4) scores for patients with Inflammatory Bowel Disease (IBD). (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2877-3)

Currently using network analysis and protein-protein interaction data to investigate and try to identify sub-networks or funtional modules that are indicative of mechanism.

Further modules and analysis will be updated in the coming months, including Cytoscape modules.

Python modules include:

`get_data.py` `get_patient.py` `get_PPI.py`


## Requirements 
* Python 3.8.x
* Genepy score matrix
* Disease Phenotype table (PatientID and Diagnosis)
* _Optional:_ LOEUF score (.xlsx) from gnomAD (https://www.nature.com/articles/s41586-020-2308-7#data-availability)
> For Python dependancies see `requirements.txt`


 
## get_data module
Functions to import Genepy matrix and generate pandas dataframes for normalised scores and subsets for Crohn's Disease and Ulcerative Colitis. Can also be used for other disease phenotypes.  
> Scores can be optionally weighted by [LOEUF](https://www.nature.com/articles/s41586-020-2308-7#data-availability) to decrease the succeptability of mutational constraint. LOEUF scores are avaliable [here](https://gnomad.broadinstitute.org/downloads) under the 'Constraint' heading.   

## get_patient module
Functions to pull down relevant individual patient data and return [StringDB](https://string-db.org/cgi/input?sessionId=bZ3itUxvQis0&input_page_active_form=single_identifier) IDs for a given list of genes from Genepy matrix using the StringDB API.  

## get_PPI module
Functions to pull down Protein-Protein Interactions (PPI's) from the StringDB using their API. 
