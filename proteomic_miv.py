import os
import urllib
import urllib.parse
import urllib.request
from bs4 import BeautifulSoup
from pandas import *

def get_uniprot (query='',query_type='PDB_ID'):
    '''
    This function will read retrieve the UniProt ID for a given 4-letter PDB file name 

    **Parameters**

    query: *string*
        The input PDB file name. query_type must be: "PDB_ID" or "ACC"
    **Returns**

        The Uniprot accession ID that corresponds to the input 4-letter PDB file name
        '''
    #code found at <a href="https://chem-workflows.com/articles/2019/10/29/retrieve-uniprot-data-using-python/">https://chem-workflows.com/articles/2019/10/29/retrieve-uniprot-data-using-python/</a>
    
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

def map_site(site_file, filename=str):
    '''
    This function will read into an csv file and convert the columns of protein. 
    Protein ID, residue start and end site number are extracted from site_file
    This function will map proteomic PSM sites in a PDB file and display in PyMOL 


    **Parameters**

    site_file: *string*
        The string of the filename for the csv file containing the proteomic site PSMs

    filename: *string*
        The string of the filename for pdb file to load into PyMOL

    **Returns**

        Text of residues with site of PSM
        PyMOL Viewer Structure with site PSM highlighted
    '''
    # reading CSV file
    site_file = "231017_PAR15map_nuclear_PSMs.csv"
    data = read_csv(site_file)

    # converting column data to list
    site_start = data['Protein.Start'].tolist()
    site_end = data['Protein.End'].tolist()
    protein_id = data['Protein.ID'].tolist()

    # create a matrix of the lists
    matrix = [protein_id, site_start, site_end]

    # access the 0 position of each of the three lists
    print("Protein ID:", matrix[0][0])
    print("Site Start:",matrix[1][0])
    print("Site End:",matrix[2][0])

    # Initialize reading into PDB file and converting into a list of lines as strings

    # To write a pml file
    bondfile=open("mapped_site.pml", "w",  encoding="utf8")
    # pylint: disable=consider-using-f-string
    bondfile.write("load {:}\n".format((filename)))
    bondfile.write("remove resn hoh\n")
    #bondfile.write("color green")

    #pdbfile=input("Enter a PDB file\n")
    rawfile=open(filename,"r",  encoding="utf8")
    pdblist=rawfile.readlines()
    rawfile.close()

    bondfile.write("color pink, resi {}-{}\n".format((matrix[1][0]), (matrix[2][0])))


### TESTING ####


#To test the get_uniprot fucntion by inputting PDB file of Solution structure of
#RRM domain in Heterogeneous nuclear ribonucleaoproteins A2/B1

#pdb_code = '1X4B'
#query_output=get_uniprot(query=pdb_code,query_type='PDB_ID')
#accession_number = query_output[1].strip().split(' ')[-1].strip(';')


#To test the download_af_pdb function by inputting an AlphaFold ID of
#Heterogeneous nuclear ribonucleoproteins A2/B1

#download_af_pdb('AF-P22626-F1')

map_site("231017_PAR15map_nuclear_PSMs.csv" , "AF-P22626-F1.pdb")