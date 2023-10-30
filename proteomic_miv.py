import os
import urllib
import urllib.parse
import urllib.request
from bs4 import BeautifulSoup

def get_uniprot (query='',query_type='PDB_ID'):
    '''
    This function will read retrieve the UniProt ID for a given 4-letter PDB file name 

    **Parameters**

    query: *string*
        The input 4-letter PDB file name
    **Returns**

        The Uniprot accession ID that corresponds to the input 4-letter PDB file name
        '''
    #code found at <a href="https://chem-workflows.com/articles/2019/10/29/retrieve-uniprot-data-using-python/">https://chem-workflows.com/articles/2019/10/29/retrieve-uniprot-data-using-python/</a>
    #query_type must be: "PDB_ID" or "ACC"
    url = 'https://www.uniprot.org/uploadlists/' #This is the webser to retrieve the Uniprot data
    params = {
    'from':query_type,
    'to':'ACC',
    'format':'txt',
    'query':query
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('ascii')
    request = urllib.request.Request(url, data)
    with urllib.request.urlopen(request) as response:
        res = response.read()
        page=BeautifulSoup(res).get_text()
        page=page.splitlines()
    return page

def uniprot_to_alphafold(accession_id):
    '''
    This function will read into the AlphaFold 'accessions_ids.csv' using the readlines command

    **Parameters**

    accession_id: *string*
        A string of the UniProt accession ID
    **Returns**

        AlphaFold ID corresponding to the the input UniProt accession ID
        '''
    accession_contents = open('accession_ids.txt', 'r').readlines()
    accession_dict = {}
    for i in accession_contents:
        accession_dict[i.split(',')[0]] = i.split(',')[3]

    accession_id = 'P14618'
    alphafold_id = None

    if accession_id in accession_dict:
        alphafold_id = accession_dict[accession_id]



def download_af_pdb(alphafold_id):
    '''
    This function will read into the AlphaFold 'accessions_ids.csv' using the readlines command

    **Parameters**

    alphafold_id: *list*
        The string of the AlphaFold ID
    **Returns**

        Download of AlphaFold PDB file
        '''
    #print("Please enter an AlphaFold ID")
    #alphafold_id = input()
    database_version = 'v4'
    model_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_id}-model_{database_version}.pdb'
    error_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_id}-predicted_aligned_error_{database_version}.json'

    os.system(f'curl {model_url} -o {alphafold_id}.pdb')
    os.system(f'curl {error_url} -o {alphafold_id}.json')

#To test the get_uniprot fucntion
'''
pdb_code = '4FXF'
query_output=get_uniprot(query=pdb_code,query_type='PDB_ID')
accession_number = query_output[1].strip().split(' ')[-1].strip(';')
'''

download_af_pdb('AF-P22626-F1')