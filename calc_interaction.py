############################################################
#####################  Core Functions  #####################
############################################################

import math

def pdb_read(pdbfile):
    '''
        This function will read into a PDB file format using the readlines command

        **Parameters**

        PDBlist: *list*
            A list of strings in PDB file format
        **Returns**

            Sorted list by acsending resiue number position
        '''
    pdbfile=input("Enter a PDB file\n")
    rawfile=open(pdbfile,"r",  encoding="utf8")
    pdblist=rawfile.readlines()
    rawfile.close()
    return pdblist

def remove(string=str):
    '''
    This function will remove space characters from a string

    **Parameters**

    string: *str*
        The given string to have removed spaces

    **Returns**

        String containing no space characters
    '''
    return string.replace(" ", "")

def only_pdb_first(string=str):
    '''
    This function will remove 'PDB_Files' from the 'filename' variable

    **Parameters**

    string: *str*
        The given string to have removed characters

    **Returns**

        String containing no 'PDB_Files'
    '''
    return string.replace("PDB_Files/", "")

def only_pdb_last(string=str):
    '''
    This function will remove '.pdb' from the 'filename' variable

    **Parameters**

    string: *str*
        The given string to have removed characters

    **Returns**

        String containing no '.pdb'
    '''
    return string.replace(".pdb", "")

def extractxyz(pdbline=str):
    '''
    This function will extract the x, y and z coordinates of an atom from the line of a PDB file

    **Parameters**

    pdbline: *str*
        

    **Returns**

        x, y and z coordinates as a list of floats
    '''
    return([float(pdbline[29:38]),float(pdbline[38:46]),float(pdbline[46:54])])

def calcdistance(list1=list,list2=list):
    '''
    This function will calculate the distance in angstroms between two atoms

    **Parameters**

    list1: *list*
        First list of floats of xyz coordniates for atoms in a PDB file
    
    list1: *list*
        Second list of floats of xyz coordniates as floats for atoms
        

    **Returns**

        Distance in angstroms between two atoms in a PDB file
    '''
    return(math.sqrt((list1[0]-list2[0])**2 + (list1[1]-list2[1])**2 + (list1[2]-list2[2])**2))

def atom_finder(pdblist=list,base=str, atom=str, mylist=list):
    '''
    This function will find a given a given atom in a nucleotide 
    base and append the whole PDB line to a given list

    **Parameters**

    base: *str*
        The nucleotide base that is to be parsed 
        (e.g. A = adenine, G = guanine, C = cytosine, T = thymine)
    
    atom: *str*
        The atom position in the nucleotide base, where the first character 
        is a letter for the elements (e.g. C=carbon, N=nitrogen, O=oxygen),
        and the second charcter is the carbon atom position in the base 
        (e.g. N7 =  the nitrogen atom on carbon 7)
    
    mylist: *list*
        Any given list to be appended with the pdbline of the atom found

    **Returns**

        Appending the atom to a given a list
    '''
    for line in pdblist:
        if line[0:4]=="ATOM":
            if line[18:20]==base:
                if line[13:16]==atom:
                    mylist.append(line)

def residue_finder(pdblist=list, resn="CYS", mylist="cys_list_unsorted"):
    '''
    This function will read PDB file and search for a given amino acid resiue, 
    and the atomic element specified. This atom will be appended to a list given

    **Parameters**

    pdbfile: *list*
        A list of strings in PDB file format
    **Returns**

        Sorted list by acsending resiue number position
    '''
    for line in pdblist:
        if (line[0:4]=="ATOM") and (line[17:20]==resn) and (line[77:78]=="S"):
            mylist.append(line)

def sorted_pdb(pdblist=list):
    '''
    This function will read a given list in PDB list, 
    and sort to return a list in acsening residue position order

    **Parameters**

    PDBlist: *list*
        A list of strings in PDB file format
    **Returns**

        Sorted list by acsending resiue number position
    '''
    return sorted(pdblist, key=lambda s: s[22:29])

def bond_distance(mylist=list, lower=1.95, upper=2.00):
    '''
    This function will calculate if a bond distance is in a given range, then print out values

    **Parameters**

    mylist: *list*
        A list of distances that correspond to potential residue-residue pairs

    lower: *float*
        The lower potential bond distance in angstroms

    uper: *float*
    The upper potential bond distance in angstroms

    **Returns*
    
        Sorted list by acsending resiue number position
    '''
    for i in mylist:
        if upper >i > lower:
            # pylint: disable=consider-using-f-string
            print("{}, {}".format(i[17:29],i))

def atom_to_atom(mylist=list, new_list=list):
    '''
    This function will go through a PDB lists of lines
    and calculatr all the potential Cysteine residue pairs

    **Parameters**

    mylist: *list*
        A list of of any amount of PDB lines. Can be specific 
        for just one amino acid residue type, or all amino acids.

    new_list: *list*
        The provided new list that PDB lines will be appended to upon the fucntion

    **Returns*
    
        Sorted list by acsending resiue number position
    '''
    # pylint: disable=consider-using-enumerate
    for i in range(len(mylist)):
        for j in range(i):
            new_list.append((mylist[i]) + (mylist[j]))

def oxygen_finder(pdblist=list, atom_name="O", new_list=list):
    '''
    This function will find oxygen atoms in a PDB list and append them to a new list

    **Parameters**

    pdblist: *list*
        A list of readlines of a PDB file

    atom_name: *str*
        The given oxygen atom to be searched (e.g. 'O')

    new_list: *list*
        The provided new list that PDB lines will be appended to upon the fucntion

    **Returns*
    
       Appended list containing all the backbone carbonyl oxygen atoms of the polypeptide
    '''
    for line in pdblist:
        if (line[0:4]=="ATOM") and (line[13:16]==atom_name):
            new_list.append(line)

def nitrogen_finder(pdblist=list, atom_name="N", new_list=list):
    '''
    This function will find oxygen atoms in a PDB list and append them to a new list

    **Parameters**

    pdblist: *list*
        A list of readlines of a PDB file

    atom_name: *str*
        The given Nitrogen atom to be searched (e.g. 'N')

    new_list: *list*
        The provided new list that PDB lines will be appended to upon the fucntion

    **Returns*
    
        Appended list containing all the backbone nitrogen atoms of the polypeptide
    '''
    for line in pdblist:
        if (line[0:4]=="ATOM") and (line[13:16]==atom_name):
            new_list.append(line)

def clean_file_path(string=str):
    '''
    This function will remove unneccessary characters in file path string

    **Parameters**

    string: *str*
        The given string to have removed spaces

    **Returns**

        String containing file path
    '''
    return string[3:-20]

def extract_pdb_path(string=str):
    '''
    This function will extract the pdb name (e.g. ifdl.pdb)
    from a given file path

    **Parameters**

    string: *str*
        The given string to have removed spaces

    **Returns**

        String containing PDB file name
    '''
    return string[-8:-4]

def extract_pdb_init(string=str):
    '''
    This function will extract the pdb name (e.g. ifdl.pdb)
    from a __init__.py to read

    **Parameters**

    string: *str*
        The given string to have removed spaces

    **Returns**

        String containing PDB file name
    '''
    return string[-4:]

############################################################
###################  Output Peptide FASTA  #################
############################################################

def output_fasta(filename=str, fasta_seq_list=list):
    '''
    This function will output a text of FASTA sequence of peptide in single amino acid code

    **Parameters**

    filename: *str*
        A string of the input PDB file

    fasta_seq_list: *list*
        An empty list to be appended to with the single amino acid code of the peptide



    **Returns*
    
        Appended list an fasta sequecne of given peptide in PDB file
    '''
    aminoacid={}
    aminoacid["ALA"]="A"
    aminoacid["CYS"]="C"
    aminoacid["ASP"]="D"
    aminoacid["GLU"]="E"
    aminoacid["PHE"]="F"
    aminoacid["GLY"]="G"
    aminoacid["HIS"]="H"
    aminoacid["ILE"]="I"
    aminoacid["LYS"]="K"
    aminoacid["LEU"]="L"
    aminoacid["MET"]="M"
    aminoacid["MSE"]="M"
    aminoacid["ASN"]="N"
    aminoacid["PRO"]="P"
    aminoacid["GLN"]="Q"
    aminoacid["ARG"]="R"
    aminoacid["SER"]="S"
    aminoacid["THR"]="T"
    aminoacid["VAL"]="V"
    aminoacid["TRP"]="W"
    aminoacid["UNK"]="X"
    aminoacid["TYR"]="Y"

    # read in a pdb file from command line:

    pdbfileraw=open(filename,"r", encoding="utf8")
    pdbfilelist=pdbfileraw.readlines()

    # make a list for storing amino acids
    aalist=[]

    for line in pdbfilelist:
        if line[0:4]=="ATOM":
            # here only look at CA lines
            if line[13:15]=="CA":
            #copy residue type into list
                aalist.append(line[17:20])

    # here output in FASTA format, with first line beginning with ">" and having info about sequence
    counter=0 # for printing out nicely
    fasta_seq_list=[]
    print(">"+filename)
    for letter in aalist:
        if letter in aminoacid:
            fasta_seq_list.append(aminoacid[letter])
            print(aminoacid[letter],end="")
        else:
            fasta_seq_list.append("X")
            print("X",end="")
    if (counter+1)%40==0:
        print()
    counter+=1

############################################################
###################  Detect Sulfide Bonds  #################
############################################################

def calc_disulfide(filename=str):
    '''
        This function will calculate any disulfide bonds in a PDB file and display in PyMOL 

        **Parameters**

        filename: *string*
            A string with the PDB file name (e.g. 1fdl.pdb)
        **Returns**

            Text of residues with connecting disulfiude bonds
            PyMOL Viewer Structure with disulfiees highlighted and bonds drawn
        '''
    # Initialize reading into PDB file and converting into a list of lines as strings

    # To write a pml file
    bondfile=open("disulfide_bonds.pml", "w",  encoding="utf8")
    # pylint: disable=consider-using-f-string
    bondfile.write("load {:}\n".format((filename)))
    bondfile.write("remove resn hoh\n")
    #bondfile.write("color green")

    #pdbfile=input("Enter a PDB file\n")
    rawfile=open(filename,"r",  encoding="utf8")
    pdblist=rawfile.readlines()
    rawfile.close()
    cys_list_unsorted=[]

    # To find all the Cysteine residues in the PDB structure
    residue_finder(pdblist,"CYS", cys_list_unsorted)

    # To sort all the Cysteine residues in acsending order
    cys_sorted=sorted_pdb(cys_list_unsorted)

    # To print out the total number of Cysteine residues in the PDB structure
    print("\nThere are",len(cys_sorted), "CYS residues")

    # A dictionary where the keys are the Cys to Cys
    # residue pairs and the values are the bond distances
    cys_dist_list=[]
    cyskeys= []
    cysvalues = []
    cysdict={}

    # To loop through the sorted Cysteine list of residues,
    # and if residues are sequential, then calculate the potential bond distances and
    # store in a new list, cys_dist_list
    # pylint: disable=consider-using-enumerate
    for i in range(len(cys_sorted)):
        for j in range(i):
            cys_dist_list.append(calcdistance(extractxyz(cys_sorted[i]), extractxyz(cys_sorted[j])))

    # To provide the the potential S-S bond interactions")
    potential_disulfide_list=[]
    atom_to_atom(cys_sorted, potential_disulfide_list)

    # To make list of Cysteine-Cysteine distances
    for i in cys_dist_list:
        cysvalues.append(i)

    # To make list of Cysteine-Cysteine combinations
    for i in potential_disulfide_list:
        cyskeys.append(i)

    # To make a dictionary matching the residue combination to the bond distance
    # pylint: disable=consider-using-enumerate
    for i in range(len(cyskeys)):
        cysdict[cyskeys[i]] = cysvalues[i]

    # To detect if a potential disulfide bond is within
    # the accpetable 2.00 angstroms plus or minus 0.05

    # pylint: disable=consider-using-dict-items
    cys_vtk = {i for i in cysdict if 2.05 >cysdict[i] > 1.95}
    true_cys_bonds_list = list(cys_vtk)
    true_cys_dist_list = list()
    for i in cysdict.values():
        if 2.05 > i> 1.95:
            true_cys_dist_list.append(i)


    # To print out the Cysteine RESN to Cysteine RESN combinations
    print("There are", len(true_cys_dist_list), "disulfide bonds\n")
    print("DISULFDE BONDS ( 2 ± 0.05 Å )")

    # To print the joined_list
    for i in true_cys_bonds_list:
        print(i[17:27], "---", i[97:107],)
    print("\nThanks for using me!")

    for line in true_cys_bonds_list:
        # pylint: disable=line-too-long
        #print("dist disulfide_bond, /{}//{}/{}`{}/{},/{}//{}/{}`{}/{}\n".format(extract_pdb_path(filename), line[21:22],line[17:20],remove(line[23:26]),line[13:16], extract_pdb_path(filename),line[102:103],remove(line[98:102]),remove(line[103:107]),line[94:96]))
        bondfile.write("dist disulfide_bond, /{}//{}/{}`{}/{},/{}//{}/{}`{}/{}\n".format(extract_pdb_path(filename), line[21:22],line[17:20],remove(line[23:26]),line[13:16], extract_pdb_path(filename),line[102:103],remove(line[98:102]),remove(line[103:107]),line[94:96]))
        bondfile.write("show sticks, /{}//{}/{}`{}\n".format(extract_pdb_path(filename), line[21:22],line[17:20],remove(line[23:26])))
        bondfile.write("show sticks, /{}//{}/{}`{}\n".format(extract_pdb_path(filename),line[102:103],remove(line[98:102]),remove(line[103:107])))
        bondfile.write("color atomic, /{}//{}/{}`{}\n".format(extract_pdb_path(filename), line[21:22],line[17:20],remove(line[23:26])))
        bondfile.write("color atomic, /{}//{}/{}`{}\n".format(extract_pdb_path(filename),line[102:103],remove(line[98:102]),remove(line[103:107])))
        bondfile.write("color atomic, disulfide")
        bondfile.write("show sticks, disulfide")
#Additional changes to alter pymol image
    bondfile.write("\nhide labels, disulfide_bond\n")
    bondfile.write("set dash_length, 0.2500\n")
    bondfile.write("set dash_gap, 0.4\n")
    bondfile.write("set dash_radius, .15\n")

############################################################
#############  WC and Non-WC Nucleic Acid Interactions  ####
############################################################

def calc_wc_nwc(filename=str):
    '''
    This function calculate all Watson-Crick 
    and Non Watson-Crick bonds in a nucleic acid PDB file

    **Parameters**

    filename: *str*
        A string of the input PDB file

    **Returns*
    
        Appended list an fasta sequecne of given peptide in PDB file
    '''
    # To write a pml file
    bondfile=open("get_bonds.pml", "w", encoding="utf8")
    # pylint: disable=consider-using-f-string
    bondfile.write("load {:}\n".format((filename)))
    bondfile.write("remove resn hoh\n")
    #bondfile.write("color green")

    rawfile = open(filename, "r", encoding="utf8")
    pdblist=rawfile.readlines()
    rawfile.close()

    def wc_dist(list1=list, list2=list):
        '''
        This function will calc distance between two lists of WC atoms and write .pml

        **Parameters**

        list1: *list*
            First list of atoms

        list2: *list*
            Second list of atoms

        **Returns*
        
            Sorted list by acsending resiue number position
        '''
        for line in list1:
            dist1=extractxyz(line)
            for i in list2:
                dist2=extractxyz(i)
                distance = calcdistance(dist1, dist2)
                if distance < 3.2:
                    # pylint: disable=line-too-long
                    bondfile.write("dist WC_hbond, /{}//{}/{}`{}/{} ,/{}//{}/{}`{}/{} , 3.2 \n".format(extract_pdb_path(filename), remove(line[21:22]),remove(line[18:20]),remove(line[23:26]),remove(line[13:16]), extract_pdb_path(filename),i[21:22],remove(i[18:20]),remove(i[23:26]),remove(i[13:15])))

    def nwc_dist(list1=list, list2=list):
        '''
        This function will calc distance between two lists of non-WC atoms and write .pml

        **Parameters**

        list1: *list*
            First list of atoms

        list2: *list*
            Second list of atoms

        **Returns*
        
            Sorted list by acsending resiue number position
        '''
        for line in list1:
            dist1=extractxyz(line)
            for i in list2:
                dist2=extractxyz(i)
                distance = calcdistance(dist1, dist2)
                if 2.5 < distance < 3.2:
                    # pylint: disable=line-too-long
                    bondfile.write("dist Non_WC_hbond, /{}//{}/{}`{}/{},/{}//{}/{}`{}/{}\n".format(extract_pdb_path(filename),line[21:22],line[19:20],remove(line[23:26]),line[13:16], extract_pdb_path(filename),i[21:22],i[19:20],i[23:26],i[13:15]))
    #lists of WC RNA bonding atoms

    #Guanine
    #pylint: disable=invalid-name
    GO6=[]
    GN1=[]
    GN2=[]

    #Cytosine
    CN4=[]
    CN3=[]
    CO2=[]

    #Adenine
    AN6=[]
    AN1=[]

    #Uracil
    UO4=[]
    UN3=[]

    #lists of WC DNA bonding atoms

    #Guanine
    #pylint: disable=invalid-name
    DGO6=[]
    DGN1=[]
    DGN2=[]

    #Cytosine
    DCN4=[]
    DCN3=[]
    DCO2=[]

    #Adenine
    DAN6=[]
    DAN1=[]

    #Uracil
    DTO4=[]
    DTN3=[]

    #lists of Non-WC bonding atoms

    #Guanine
    GN3=[]
    GN9=[]
    GN7=[]

    #Adenine
    AN7=[]
    AN9=[]
    AN3=[]

    #Uracil
    UO2=[]

    #WC RNA atoms
    (atom_finder(pdblist,"G","O6 ",GO6))
    (atom_finder(pdblist,"G","N1 ",GN1))
    (atom_finder(pdblist,"G","N2 ",GN2))

    (atom_finder(pdblist,"C","N4 ",CN4))
    (atom_finder(pdblist,"C","N3 ",CN3))
    (atom_finder(pdblist,"C","O2 ",CO2))

    (atom_finder(pdblist,"A","N6 ",AN6))
    (atom_finder(pdblist,"A","N1 ",AN1))

    (atom_finder(pdblist,"U","O4 ",UO4))
    (atom_finder(pdblist,"U","N3 ",UN3))

    #WC DNA atoms
    (atom_finder(pdblist,"DG","O6 ",DGO6))
    (atom_finder(pdblist,"DG","N1 ",DGN1))
    (atom_finder(pdblist,"DG","N2 ",DGN2))

    (atom_finder(pdblist,"DC","N4 ",DCN4))
    (atom_finder(pdblist,"DC","N3 ",DCN3))
    (atom_finder(pdblist,"DC","O2 ",DCO2))

    (atom_finder(pdblist,"DA","N6 ",DAN6))
    (atom_finder(pdblist,"DA","N1 ",DAN1))

    (atom_finder(pdblist,"DT","O4 ",DTO4))
    (atom_finder(pdblist,"DT","N3 ",DTN3))



    #Non-WC RNA atoms
    (atom_finder(pdblist,"G","N3 ",GN3))
    (atom_finder(pdblist,"G","N9 ",GN9))
    (atom_finder(pdblist,"G","N7 ",GN7))


    (atom_finder(pdblist,"A","N7 ",AN7))
    (atom_finder(pdblist,"A","N9 ",AN9))
    (atom_finder(pdblist,"A","N3 ",AN3))

    (atom_finder(pdblist,"U","O2 ",UO2))

    # WC RNA write .pml
    wc_dist(GO6, CN4)
    wc_dist(GN1, CN3)
    wc_dist(GN2, CO2)

    wc_dist(AN6, UO4)
    wc_dist(AN1, UN3)

    # WC DNA write .pml
    wc_dist(DGO6, DCN4)
    wc_dist(DGN1, DCN3)
    wc_dist(DGN2, DCO2)

    wc_dist(DAN6, DTO4)
    wc_dist(DAN1, DTN3)
    bondfile.write("set dash_color, yellow, WC_hbond")

    #Non-WC distances
    #Series of all the atoms tested for hydrogen bond (these are all non-WC)
    print(DAN6)
    print(AN6)
    if not DAN6:
        print("Hellow world")
        #Hoogsteen
        nwc_dist(AN7,UN3)
        nwc_dist(AN6,UO4)
        nwc_dist(GN7,CN3)
        #GO6
        nwc_dist(GO6,CO2)
        nwc_dist(GO6,AN6)
        nwc_dist(GO6,AN1)
        nwc_dist(GO6,AN7)
        nwc_dist(GO6,UN3)
        #GN1
        nwc_dist(GN1,GN3)
        nwc_dist(GN1,GN7)
        nwc_dist(GN1,AN1)
        nwc_dist(GN1,AN3)
        nwc_dist(GN1,AN7)
        nwc_dist(GN1,UO4)
        nwc_dist(GN1,UN3)
        #GN2
        nwc_dist(GN2,GN7)
        nwc_dist(GN2,AN6)
        nwc_dist(GN2,AN1)
        nwc_dist(GN2,AN3)
        nwc_dist(GN2,AN7)
        nwc_dist(GN2,UO4)
        nwc_dist(GN2,UN3)
        #GN3
        nwc_dist(GN3,GO6)
        nwc_dist(GN3,CN3)
        nwc_dist(GN3,CO2)
        nwc_dist(GN3,AN6)
        nwc_dist(GN3,AN1)
        nwc_dist(GN3,AN3)
        nwc_dist(GN3,UO4)
        nwc_dist(GN3,UN3)
        #GN7
        nwc_dist(GN7,GN1)
        nwc_dist(GN7,GN2)
        nwc_dist(GN7,GN7)
        nwc_dist(GN7,CN3)
        nwc_dist(GN7,CO2)
        nwc_dist(GN7,AN6)
        nwc_dist(GN7,AN1)
        nwc_dist(GN7,AN3)
        nwc_dist(GN7,AN7)
        nwc_dist(GN7,UO4)
        nwc_dist(GN7,UN3)
        #AN6
        nwc_dist(AN6,GO6)
        nwc_dist(AN6,GN2)
        nwc_dist(AN6,GN3)
        nwc_dist(AN6,GN7)
        nwc_dist(AN6,CN3)
        nwc_dist(AN6,CO2)
        nwc_dist(AN6,UN3)
        #AN1
        nwc_dist(AN1,GO6)
        nwc_dist(AN1,GN1)
        nwc_dist(AN1,GN2)
        nwc_dist(AN1,GN3)
        nwc_dist(AN1,GN7)
        nwc_dist(AN1,CN3)
        nwc_dist(AN1,CO2)
        nwc_dist(AN1,AN6)
        nwc_dist(AN1,AN1)
        nwc_dist(AN1,AN3)
        nwc_dist(AN1,AN7)
        nwc_dist(AN1,UO4)
        #AN7
        nwc_dist(AN7,GO6)
        nwc_dist(AN7,GN1)
        nwc_dist(AN7,GN2)
        nwc_dist(AN7,GN7)
        nwc_dist(AN7,CN3)
        nwc_dist(AN7,CO2)
        nwc_dist(AN7,AN1)
        nwc_dist(AN7,AN3)
        nwc_dist(AN7,UO4)
        nwc_dist(AN7,UN3)
        #AN3
        nwc_dist(AN3,GO6)
        nwc_dist(AN3,GN1)
        nwc_dist(AN3,GN2)
        nwc_dist(AN3,GN3)
        nwc_dist(AN3,GN7)
        nwc_dist(AN3,CN3)
        nwc_dist(AN3,CO2)
        nwc_dist(AN3,AN6)
        nwc_dist(AN3,AN7)
        nwc_dist(AN3,AN9)
        nwc_dist(AN3,UO4)
        nwc_dist(AN3,UN3)
        #UO4
        nwc_dist(UO4,GN1)
        nwc_dist(UO4,GN2)
        nwc_dist(UO4,GN3)
        nwc_dist(UO4,GN7)
        nwc_dist(UO4,CN3)
        nwc_dist(UO4,CO2)
        nwc_dist(UO4,AN1)
        nwc_dist(UO4,AN3)
        nwc_dist(UO4,AN7)
        #UN3
        nwc_dist(UN3,GO6)
        nwc_dist(UN3,GN1)
        nwc_dist(UN3,GN2)
        nwc_dist(UN3,GN3)
        nwc_dist(UN3,GN7)
        nwc_dist(UN3,CN3)
        nwc_dist(UN3,CO2)
        nwc_dist(UN3,AN6)
        nwc_dist(UN3,AN3)
        nwc_dist(UN3,AN7)
        nwc_dist(UN3,UO4)
        nwc_dist(UO4,UO2)
        #CN4
        nwc_dist(CN4,GN3)
        nwc_dist(CN4,GN7)
        nwc_dist(CN4,AN1)
        nwc_dist(CN4,AN3)
        nwc_dist(CN4,AN7)
        nwc_dist(CN4,UO4)
        nwc_dist(CN4,UN3)
        #CN3
        nwc_dist(CN3,GN3)
        nwc_dist(CN3,GN7)
        nwc_dist(CN3,AN6)
        nwc_dist(CN3,AN1)
        nwc_dist(CN3,AN3)
        nwc_dist(CN3,AN7)
        nwc_dist(CN3,UO4)
        nwc_dist(CN3,UN3)
        #CO2
        nwc_dist(CO2,GO6)
        nwc_dist(CO2,GN3)
        nwc_dist(CO2,GN7)
        nwc_dist(CO2,AN6)
        nwc_dist(CO2,AN1)
        nwc_dist(CO2,AN3)
        nwc_dist(CO2,AN7)
        nwc_dist(CO2,UO4)
        nwc_dist(CO2,UN3)

    #Additional changes to alter pymol image
    bondfile.write("\n")
    bondfile.write("set dash_color, red, Non_WC_hbond")
    bondfile.write("\n")
    bondfile.write("set cartoon_ring_mode,3")
    bondfile.write("\n")
    bondfile.write("hide labels, Non_WC_hbond")
    bondfile.write("\n")
    bondfile.write("hide labels, WC_hbond")
    bondfile.write("\n")
    bondfile.write("show sticks, sidechain extend 1")
    bondfile.write("\n")
    bondfile.write("color grey, sidechain")
    bondfile.write("\n")
    bondfile.write("color atomic, (not elem C)")
    bondfile.write("\n")
    bondfile.write("set dash_length, 0.2500")
    bondfile.write("\n")
    bondfile.write("set dash_gap, 0.4")
    bondfile.write("\n")
    bondfile.write("set dash_radius, .15")

############################################################
###################  Detect Alpha Helice  ##################
############################################################

def alpha_helice(filename=str):
    '''
        This function will detect any potential alpha helical seconday structure in polypeptides

        **Parameters**

        filename: *string*
            A string with the PDB file name (e.g. 1fdl.pdb)
        **Returns**

            Text of residues with alpha helical structure
            PyMOL Viewer Structure with alpha helices highlighted and bonds drawn
        '''
    # To write a pml file
    bondfile=open("helix_bonds.pml", "w", encoding="utf8")
    # pylint: disable=consider-using-f-string
    bondfile.write("load {:}\n".format((filename)))
    bondfile.write("remove resn hoh\n")
    bondfile.write("color white\n")

    aminoacid={}
    aminoacid["ALA"]="A"
    aminoacid["CYS"]="C"
    aminoacid["ASP"]="D"
    aminoacid["GLU"]="E"
    aminoacid["PHE"]="F"
    aminoacid["GLY"]="G"
    aminoacid["HIS"]="H"
    aminoacid["ILE"]="I"
    aminoacid["LYS"]="K"
    aminoacid["LEU"]="L"
    aminoacid["MET"]="M"
    aminoacid["MSE"]="M"
    aminoacid["ASN"]="N"
    aminoacid["PRO"]="P"
    aminoacid["GLN"]="Q"
    aminoacid["ARG"]="R"
    aminoacid["SER"]="S"
    aminoacid["THR"]="T"
    aminoacid["VAL"]="V"
    aminoacid["TRP"]="W"
    aminoacid["UNK"]="X"
    aminoacid["TYR"]="Y"

    #to put all the carbonly oxygen ATOM lines in a list
    oxy_atom_list=[]
    nitro_atom_list=[]

    # To read into PDB file
    print("Alpha-helical structure detector\n")
    rawfile=open(filename,"r", encoding="utf8")
    pdblist=rawfile.readlines()
    rawfile.close()


    oxygen_finder(pdblist,"O  ", oxy_atom_list)
    #print(oxy_atom_list)

    nitrogen_finder(pdblist,"N  ", nitro_atom_list)
    #print(nitro_atom_list)

    # To remove the extraneous last 4 nitrogens that cannot participate in alpha helix
    try:
        nitro_atom_list.pop(0)
        nitro_atom_list.pop(0)
        nitro_atom_list.pop(0)
        nitro_atom_list.pop(0)
    except IndexError:
        print("Please enter a valid protein PDB file")


    dist_list=[]
    # pylint: disable=consider-using-enumerate
    for i in range(len(nitro_atom_list)):
        dist_list.append(calcdistance(extractxyz(oxy_atom_list[i]), extractxyz(nitro_atom_list[i])))
    aalist=[]
    for line in oxy_atom_list:
        aalist.append(line[17:20])

    # To output the protein sequence in FASTA format
    single_aa_list=[]

    #Make a list of the single letter amino acid FASTA sequence
    for letter in aalist:
        if letter in aminoacid:
            single_aa_list.append(aminoacid[letter])
        else:
            single_aa_list.append("x")

    #Make a list of the "-" and "H" for the single amino acid FASTA seqeuence
    h_bond_list=[]
    for i in range(len(dist_list)):
        if 0.5 < dist_list[i] < 4.7:
            h_bond_list.append("H")
        else:
            h_bond_list.append("-")
    helix_line_list_n=[]
    helix_line_list_o=[]
    for i in range(len(nitro_atom_list)):
        if 0.5 < (calcdistance(extractxyz(oxy_atom_list[i]), extractxyz(nitro_atom_list[i]))) < 4.7:
            helix_line_list_n.append(nitro_atom_list[i])
            helix_line_list_o.append(oxy_atom_list[i])

    print(helix_line_list_n)
    print("\n\n\n")
    print(helix_line_list_o)

    # To account for the 4 i+4 Hydorgen bonding between carbonyl oxygen and amine
    h_bond_list.append("-")
    h_bond_list.append("-")
    h_bond_list.append("-")
    h_bond_list.append("-")

    print("\n\n'H' = alpha helical structure")
    print("'-' = non-alpha helical structure")
    print("\n>"+filename)

    #Here is a method to print the fasta sequence and alpha helical structure
    #in an alinged manner, coded up to 531 characters, so this method works for aa < 531, and
    #can be coded for higher MW proteins if necessary

    for line in helix_line_list_n:
        # pylint: disable=line-too-long
        #bondfile.write("show sticks, /{}//{}/{}`{}\n".format(extract_pdb_path(filename),line[102:103],remove(line[98:102]),remove(line[103:107])))
        bondfile.write("color pink, /{}//{}/{}`{}\n".format(extract_pdb_path(filename), line[21:22],line[17:20],remove(line[23:26])))


    for line in helix_line_list_o:
        # pylint: disable=line-too-long
        #bondfile.write("show sticks, /{}//{}/{}`{}\n".format(extract_pdb_path(filename),line[102:103],remove(line[98:102]),remove(line[103:107])))
        bondfile.write("color pink, /{}//{}/{}`{}\n".format(extract_pdb_path(filename), line[21:22],line[17:20],remove(line[23:26])))
#Additional changes to alter pymol image

    print(*single_aa_list[0:40], sep = "")
    print(*h_bond_list[0:40], sep = "")

    print(*single_aa_list[41:81], sep = "")
    print(*h_bond_list[41:81], sep = "")

    print(*single_aa_list[82:122], sep = "")
    print(*h_bond_list[82:122], sep = "")

    print(*single_aa_list[163:203], sep = "")
    print(*h_bond_list[163:203], sep = "")

    print(*single_aa_list[204:244], sep = "")
    print(*h_bond_list[204:244], sep = "")

    print(*single_aa_list[245:285], sep = "")
    print(*h_bond_list[245:285], sep = "")

    print(*single_aa_list[286:326], sep = "")
    print(*h_bond_list[286:326], sep = "")

    print(*single_aa_list[327:367], sep = "")
    print(*h_bond_list[327:367], sep = "")

    print(*single_aa_list[368:408], sep = "")
    print(*h_bond_list[368:408], sep = "")

    print(*single_aa_list[409:449], sep = "")
    print(*h_bond_list[409:449], sep = "")

    print(*single_aa_list[450:490], sep = "")
    print(*h_bond_list[450:490], sep = "")

    print(*single_aa_list[491:531], sep = "")
    print(*h_bond_list[491:531], sep = "")

    print(*single_aa_list[532:572], sep = "")
    print(*h_bond_list[532:572], sep = "")

############################################################
#####################  End to End Distance  ################
############################################################

def end_to_end_dist(filename=str):
    '''
        This function will calculate the molecular weight of a peptide

        **Parameters**

        filename: *string*
            A string with the PDB file name (e.g. 1fdl.pdb)
        **Returns**

            Text of molecular weight of a peptide 
        '''

    # To write a pml file
    bondfile=open("end_to_end.pml", "w", encoding="utf8")
    # pylint: disable=consider-using-f-string
    bondfile.write("load {:}\n".format((filename)))
    bondfile.write("remove resn hoh\n")

    while True:
    #will read into file if ".pdb"
        try:
            rawfile=open(filename,"r", encoding="utf8")
            pdblist=rawfile.readlines()
            rawfile.close()
            rawfile=open(filename,"r", encoding="utf8")
            pdblist=rawfile.readlines()
            rawfile.close()

            #Generate a file name to store CA atoms
            outputfile= "outputMWFile"
            ca_atoms=open(outputfile,"w", encoding="utf8")

            def strip_ca_atoms(pdb):
                ### goal of function is to take PDB list and return only CA atoms
                returnedlist=[]
                for line in pdb:
                    if line[0:4]=="ATOM" and line[13:15]=="CA":
                        ca_atoms.write(line)
                        returnedlist.append(line)
                return returnedlist
            # pylint: disable= unused-variable
            ca_only=strip_ca_atoms(pdblist)

            ca_atoms.close()
            #extract the first and last atoms from CA atoms file
            with open(outputfile, "r",  encoding="utf8") as file:
                first_line = file.readline()
                for last_line in file:
                    pass
            # pylint: disable=undefined-loop-variable
            ca_only = [first_line, last_line]

            # Writing to end_to_end.pml
            #pylint: disable=line-too-long
            bondfile.write("show sticks, /{}//{}/{}`{}\n".format(extract_pdb_path(filename),last_line[21:22],last_line[17:20],last_line[23:26]))
            bondfile.write("dist end_to_end, /{}//{}/{}`{}/{},/{}//{}/{}`{}/{}\n".format(extract_pdb_path(filename), first_line[21:22],first_line[17:20],remove(first_line[23:26]),first_line[13:16], extract_pdb_path(filename),last_line[21:22],last_line[17:20],last_line[23:26],last_line[13:15]))
            bondfile.write("show sticks, /{}//{}/{}`{}\n".format(extract_pdb_path(filename), first_line[21:22],first_line[17:20],remove(first_line[23:26])))
            bondfile.write("show sticks, /{}//{}/{}`{}\n".format(extract_pdb_path(filename),last_line[21:22],last_line[17:20],last_line[23:26]))
            bondfile.write("set dash_color, pink, end_to_end\n")
            #bondfile.write("hide labels, Non_WC_hbond")
            #bondfile.write("color atomic, (not elem C), /{}//{}/{}`{}/{}\n".format(extract_pdb_path(filename), first_line[21:22],first_line[19:20],remove(first_line[23:26]),first_line[13:16]))
            #bondfile.write("color atomic, (not elem C), /{}//{}/{}`{}/{}\n".format(extract_pdb_path(filename),last_line[21:22],last_line[19:20],last_line[23:26],last_line[13:15]))
            bondfile.write("set dash_length, 0.2500\n")
            bondfile.write("set dash_gap, 0.4\n")
            bondfile.write("set dash_radius, .15\n")

            # Print Statments
            print("\nFirst Residue:",first_line[17:20])
            print("Last Residue:", last_line[17:20])
            # pylint: disable=consider-using-f-string
            #pylint: disable=line-too-long
            print("Distance between Cα atoms of first and last residue: {:.2f} Å".format((calcdistance(extractxyz(first_line), extractxyz(last_line)))))
            return None
        except TypeError:
            print("Enter a valid PDB File")
            break

############################################################
###################  Calculate Peptide MW  #################
############################################################

def calc_peptide_mw(filename):
    '''
        This function will calculate the molecular weight of a peptide

        **Parameters**

        filename: *string*
            A string with the PDB file name (e.g. 1fdl.pdb)
        **Returns**

            Text of molecular weight of a peptide 
        '''
# To read into a PDB file and then outputs the protein sequence in FASTA format

    aminoacid={}
    aminoacid["ALA"]="A"
    aminoacid["CYS"]="C"
    aminoacid["ASP"]="D"
    aminoacid["GLU"]="E"
    aminoacid["PHE"]="F"
    aminoacid["GLY"]="G"
    aminoacid["HIS"]="H"
    aminoacid["ILE"]="I"
    aminoacid["LYS"]="K"
    aminoacid["LEU"]="L"
    aminoacid["MET"]="M"
    aminoacid["MSE"]="M"
    aminoacid["ASN"]="N"
    aminoacid["PRO"]="P"
    aminoacid["GLN"]="Q"
    aminoacid["ARG"]="R"
    aminoacid["SER"]="S"
    aminoacid["THR"]="T"
    aminoacid["VAL"]="V"
    aminoacid["TRP"]="W"
    aminoacid["UNK"]="X"
    aminoacid["TYR"]="Y"

    # read in a pdb file from command line:

    pdbfileraw=open(filename,"r",  encoding="utf8")
    pdbfilelist=pdbfileraw.readlines()

    # make a list for storing amino acids
    aalist=[]

    for line in pdbfilelist:
        if line[0:4]=="ATOM":
            # here only look at CA lines
            if line[13:15]=="CA":
            #copy residue type into list
                aalist.append(line[17:20])

    # here output in FASTA format, with first line beginning with ">" and having info about sequence
    counter=0 # for printing out nicely
    fasta_seq_list=[]
    print(">"+filename)
    for letter in aalist:
        if letter in aminoacid:
            fasta_seq_list.append(aminoacid[letter])
        else:
            fasta_seq_list.append("X")
    if (counter+1)%40==0:
        print()
    counter+=1

    #Amino Acid Molecular Weight dictionary
    aa_mw={}
    aa_mw["A"]=89.09
    aa_mw["C"]=121.16
    aa_mw["D"]=133.10
    aa_mw["E"]=147.13
    aa_mw["F"]=165.19
    aa_mw["G"]=75.07
    aa_mw["H"]=155.16
    aa_mw["I"]=131.18
    aa_mw["K"]=146.19
    aa_mw["L"]=131.18
    aa_mw["M"]=149.21
    aa_mw["N"]=132.12
    aa_mw["P"]=115.13
    aa_mw["Q"]=146.15
    aa_mw["R"]=174.20
    aa_mw["S"]=105.09
    aa_mw["T"]=119.12
    aa_mw["V"]=117.15
    aa_mw["W"]=204.23
    aa_mw["Y"]=181.19
    aa_mw["["]=0
    aa_mw["]"]=0
    aa_mw["\n"]=0
    aa_mw[","]=0
    aa_mw[" "]=0
    aa_mw["'"]=0
    aa_mw["_"]=0
    aa_mw["X"]=0

    #water molecular weight
    aa_mw["HOH"]=18.02

    # To print the amino acid sequence from fasta file

    total=0
    for char in fasta_seq_list:
        letter_to_number = aa_mw[char]
        total += letter_to_number

    loss_of_water = (int(len(fasta_seq_list))-1) * 18.02
    peptide_mass=total - loss_of_water
    # pylint: disable=consider-using-f-string
    print("Peptide Mass: {:.2f} Daltons".format(peptide_mass))
