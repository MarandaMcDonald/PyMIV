# To initialize the plugin to be installed by PyMOL
from __future__ import absolute_import
from __future__ import print_function

# To provide an entry point to PyMOL's API
import os
from pymol import cmd
from pymol.Qt import QtWidgets
from pymol.Qt.utils import loadUi
from pymol.Qt.utils import getSaveFileNameWithExt

# To load the UI file into our dialog
from pymol.Qt.utils import loadUi

#This is a series of functions for mathematical calculations to read into PDB files and output a .pml text to import into the PyMOL


import math

############################################################
#######################Core Functions#######################
############################################################

def PDB_read(PDBfile):
    '''
        This function will read into a PDB file format using the readlines command

        **Parameters**

        PDBlist: *list*
            A list of strings in PDB file format
        **Returns**

            Sorted list by acsending resiue number position
        '''
    PDBfile=input("Enter a PDB file\n")
    rawfile=open(PDBfile,"r")
    pdblist=rawfile.readlines()
    rawfile.close()
    return None

def remove(string):
    '''
    This function will remove space characters from a string

    **Parameters**

    string: *str*
        The given string to have removed spaces

    **Returns**

        String containing no space characters
    '''
    return string.replace(" ", "")
    
def extractxyz(PDBline):
    '''
    This function will extract the x, y and z coordinates of an atom from the line of a PDB file

    **Parameters**

    PDBline: *str*
        

    **Returns**

        x, y and z coordinates as a list of floats
    '''
    return([float(PDBline[29:38]),float(PDBline[38:46]),float(PDBline[46:54])])

def calcdistance(list1,list2):
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

def atom_finder(pdblist,base, atom, mylist):
    '''
    This function will find a given a given atom in a nucleotide base and append the whole PDB line to a given list

    **Parameters**

    base: *str*
        The nucleotide base that is to be parsed (e.g. A = adenine, G = guanine, C = cytosine, T = thymine)
    
    atom: *str*
        The atom position in the nucleotide base, where the first characte is a letter for the elements (e.g. C=carbon, N=nitrogen, O=oxygen),
        and the second charcter is the carbon atom position in the base (e.g. N7 =  the nitrogen atom on carbon 7)
    
    mylist: *list*
        Any given list to be appended with the PDBline of the atom found

    **Returns**

        Appending the atom to a given a list
    '''
    for line in pdblist:
        if (line[0:4]=="ATOM"):
            if (line[19:20]==base):
                if (line[13:16]==atom):
                    mylist.append(line)

def residue_finder(pdblist, resn="CYS", mylist="CYSlistUnsorted"):
    '''
    This function will read PDB file and search for a given amino acid resiue, and the atomic element specified. This atom will be appended to a list given

    **Parameters**

    PDBfile: *list*
        A list of strings in PDB file format
    **Returns**

        Sorted list by acsending resiue number position
    '''
    for line in pdblist:
        if (line[0:4]=="ATOM") and (line[17:20]==resn) and (line[77:78]=="S"):
            mylist.append(line)

def sorted_pdb(pdblist):
    '''
    This function will read a given list in PDB list, and sort to return a list in acsening residue position order

    **Parameters**

    PDBlist: *list*
        A list of strings in PDB file format
    **Returns**

        Sorted list by acsending resiue number position
    '''
    return(sorted(pdblist, key=lambda s: s[22:29]))

def bond_distance(mylist, lower=1.95, upper=2.00):
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
            print("{}, {}".format(i[17:29],i))

def atom_to_atom(mylist, newList):
    '''
    This function will go through a PDB lists of lines and calculatr all the potential Cysteine residue pairs

    **Parameters**

    mylist: *list*
        A list of of any amount of PDB lines. Can be specific for just one amino acid residue type, or all amino acids.

    newList: *list*
        The provided new list that PDB lines will be appended to upon the fucntion

    **Returns*
    
        Sorted list by acsending resiue number position
    '''
    for i in range(len(mylist)):
        for z in range(i):
            newList.append((mylist[i]) + (mylist[z]))

############################################################
#####################Detect Sulfide Bonds###################
############################################################

def calc_disulfide(filename):
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
    
    #PDBfile=input("Enter a PDB file\n")
    rawfile=open(filename,"r")
    pdblist=rawfile.readlines()
    rawfile.close()
    CYSlistUnsorted=[]
    
    # To find all the Cysteine residues in the PDB structure
    residue_finder(pdblist,"CYS", CYSlistUnsorted)

    # To sort all the Cysteine residues in acsending order
    CYSlist=sorted_pdb(CYSlistUnsorted)

    # To print out the total number of Cysteine residues in the PDB structure
    print("\nThere are",len(CYSlist), "CYS residues")

    # A dictionary where the keys are the Cys to Cys residue pairs and the values are the bond distances
    CYSdistanceslist=[]
    cyskeys= []
    cysvalues = []
    cysdict={}

    # To loop through the sorted Cysteine list of residues, and if residues are sequential, then calculate the potential bond distances and
    # store in a new list, CYSdistanceslist
    for i in range(len(CYSlist)):
        for z in range(i):
            CYSdistanceslist.append(calcdistance(extractxyz(CYSlist[i]), extractxyz(CYSlist[z])))

    # To provide the the potential S-S bond interactions")
    StoSlist=[]
    atom_to_atom(CYSlist, StoSlist)

    # To make list of Cysteine-Cysteine distances
    for i in CYSdistanceslist:
        cysvalues.append(i)

    # To make list of Cysteine-Cysteine combinations
    for i in StoSlist:
        cyskeys.append(i)

    # To make a dictionary matching the residue combination to the bond distance
    for i in range(len(cyskeys)):
        cysdict[cyskeys[i]] = cysvalues[i]

    # To detect if a potential disulfide bond is within the accpetable 2.00 angstroms plus or minus 0.05
    cysValueToKey = {i for i in cysdict if 2.05 >cysdict[i] > 1.95}
    trueCysBondlist = list(cysValueToKey)
    trueDistancelist = list()
    for i in cysdict.values():
            if 2.05 > i> 1.95:
                trueDistancelist.append(i)

    # To print out the Cysteine RESN to Cysteine RESN combinations
    joinedList = "\n".join("{} {}".format(x, y) for x, y in zip(trueCysBondlist, "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"))
    print("There are", len(trueDistancelist), "disulfide bonds")
    print()
    print("DISULFDE BONDS ( 2 ± 0.05 Å )")

    # To print the joinedList
    for i in trueCysBondlist:
        print(i[17:27], "---", i[97:107],)
    print()
    print("Thanks for using me!")

'''
############################################################
###############WC and Non-WC Nucleic Acid Interactions######
############################################################

#to write a pml file
bondfile=open("get_bonds.pml", "w")
bondfile.write("load 1z43.pdb\n")
bondfile.write("remove resn hoh\n")
#bondfile.write("color green")

#x=input("Enter a PDB file\n")
x='1z43.pdb'
rawfile=open(x,"r")
pdblist=rawfile.readlines()
rawfile.close()

#lists of WC bonding atoms

#Guanine
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

dist1=[]
#WC atoms
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

#Non-WC atoms
(atom_finder(pdblist,"G","N3 ",GN3))
(atom_finder(pdblist,"G","N9 ",GN9))
(atom_finder(pdblist,"G","N7 ",GN7))


(atom_finder(pdblist,"A","N7 ",AN7))
(atom_finder(pdblist,"A","N9 ",AN9))
(atom_finder(pdblist,"A","N3 ",AN3))

(atom_finder(pdblist,"U","O2 ",UO2))


#function to write to a pml file to draw the WC hydrogen bond distance for each WC hydrogen bonding atom pair
def WC_distance_pymol(list1, list2):
    for line in list1:
        dist1=extractxyz(line)
        for i in list2:
            dist2=extractxyz(i)
            distance = (calcdistance(dist1, dist2))
            if distance < 3.2:
                        bondfile.write("dist WC_hbond, /1Z43//{}/{}`{}/{},/1Z43//{}/{}`{}/{}\n".format(line[21:22],line[19:20],remove(line[23:26]),line[13:16],i[21:22],i[19:20],i[23:26],i[13:15]))
                        #bondfile.write("show sticks, /1Z43//{}/{}`{}\n".format(line[21:22],line[19:20],remove(line[23:26]),line[13:16]))
                        #bondfile.write("show sticks, /1Z43//{}/{}`{}\n".format(i[21:22],i[19:20],i[23:26],i[13:15]))

#write to pml file of 
WC_distance_pymol(GO6, CN4)
WC_distance_pymol(GN1, CN3)
WC_distance_pymol(GN2, CO2)

WC_distance_pymol(AN6, UO4)
WC_distance_pymol(AN1, UN3)

bondfile.write("set dash_color, yellow, WC_hbond")

#Non-WC distances
#function to write to a pml file to draw the WC hydrogen bond distance for each WC hydrogen bonding atom pair

def Non_WC_distance_pymol(list1, list2):
    for line in list1:
        dist1=extractxyz(line)
        for i in list2:
            dist2=extractxyz(i)
            distance = (calcdistance(dist1, dist2))
            if 2.5 < distance < 3.2:
                        bondfile.write("dist Non_WC_hbond, /1Z43//{}/{}`{}/{},/1Z43//{}/{}`{}/{}\n".format(line[21:22],line[19:20],remove(line[23:26]),line[13:16],i[21:22],i[19:20],i[23:26],i[13:15]))
                        #bondfile.write("show sticks, /1Z43//{}/{}`{}\n".format(line[21:22],line[19:20],remove(line[23:26]),line[13:16]))
                        #bondfile.write("show sticks, /1Z43//{}/{}`{}\n".format(i[21:22],i[19:20],i[23:26],i[13:15]))

atomList=[GO6,GN1,GN2,GN3,GN7,GN9,CN4,CN3,CO2,AN6,AN1,AN7,AN9,AN3,UO4,UN3,]
#WC BONDING ATOMS LIST
Gwc=[GO6,GN1,GN2]
Cwc=[CN4,CN3,CO2]
Awc=[AN6,AN1]
Uwc=[UO4,UN3]

#HOOGSTEIN/NONWC BONDING ATOMS LIST
Ghoog=[GN3,GN9,GN7]
Ahoog=[AN7,AN3,AN9]
'''
'''
#Series of all the atoms tested for hydrogen bond (these are all non-WC)
#Manually exmained the pymol structure to elimate which distances are not true hydrogen bonds
#Hoogsteen
Non_WC_distance_pymol(AN7,UN3)
Non_WC_distance_pymol(AN6,UO4)

Non_WC_distance_pymol(GN7,CN3)

#GO6
Non_WC_distance_pymol(GO6,CO2)
Non_WC_distance_pymol(GO6,AN6)
Non_WC_distance_pymol(GO6,AN1)
Non_WC_distance_pymol(GO6,AN7)

Non_WC_distance_pymol(GO6,UN3)

#GN1
Non_WC_distance_pymol(GN1,GN3)
Non_WC_distance_pymol(GN1,GN7)


Non_WC_distance_pymol(GN1,AN1)
Non_WC_distance_pymol(GN1,AN3)
Non_WC_distance_pymol(GN1,AN7)

Non_WC_distance_pymol(GN1,UO4)
Non_WC_distance_pymol(GN1,UN3)


#GN2
Non_WC_distance_pymol(GN2,GN7)

Non_WC_distance_pymol(GN2,AN6)
Non_WC_distance_pymol(GN2,AN1)
Non_WC_distance_pymol(GN2,AN3)
Non_WC_distance_pymol(GN2,AN7)

Non_WC_distance_pymol(GN2,UO4)
Non_WC_distance_pymol(GN2,UN3)


#GN3
Non_WC_distance_pymol(GN3,GO6)

Non_WC_distance_pymol(GN3,CN3)
Non_WC_distance_pymol(GN3,CO2)

Non_WC_distance_pymol(GN3,AN6)
Non_WC_distance_pymol(GN3,AN1)
Non_WC_distance_pymol(GN3,AN3)

Non_WC_distance_pymol(GN3,UO4)
Non_WC_distance_pymol(GN3,UN3)

#GN7
Non_WC_distance_pymol(GN7,GN1)
Non_WC_distance_pymol(GN7,GN2)
Non_WC_distance_pymol(GN7,GN7)

Non_WC_distance_pymol(GN7,CN3)
Non_WC_distance_pymol(GN7,CO2)

Non_WC_distance_pymol(GN7,AN6)
Non_WC_distance_pymol(GN7,AN1)
Non_WC_distance_pymol(GN7,AN3)
Non_WC_distance_pymol(GN7,AN7)

Non_WC_distance_pymol(GN7,UO4)
Non_WC_distance_pymol(GN7,UN3)


#AN6
Non_WC_distance_pymol(AN6,GO6)
Non_WC_distance_pymol(AN6,GN2)
Non_WC_distance_pymol(AN6,GN3)
Non_WC_distance_pymol(AN6,GN7)

Non_WC_distance_pymol(AN6,CN3)
Non_WC_distance_pymol(AN6,CO2)

Non_WC_distance_pymol(AN6,UN3)


#AN1
Non_WC_distance_pymol(AN1,GO6)
Non_WC_distance_pymol(AN1,GN1)
Non_WC_distance_pymol(AN1,GN2)
Non_WC_distance_pymol(AN1,GN3)
Non_WC_distance_pymol(AN1,GN7)

Non_WC_distance_pymol(AN1,CN3)
Non_WC_distance_pymol(AN1,CO2)

Non_WC_distance_pymol(AN1,AN6)
Non_WC_distance_pymol(AN1,AN1)
Non_WC_distance_pymol(AN1,AN3)
Non_WC_distance_pymol(AN1,AN7)

Non_WC_distance_pymol(AN1,UO4)

#AN7
Non_WC_distance_pymol(AN7,GO6)
Non_WC_distance_pymol(AN7,GN1)
Non_WC_distance_pymol(AN7,GN2)
Non_WC_distance_pymol(AN7,GN7)

Non_WC_distance_pymol(AN7,CN3)
Non_WC_distance_pymol(AN7,CO2)

Non_WC_distance_pymol(AN7,AN1)
Non_WC_distance_pymol(AN7,AN3)

Non_WC_distance_pymol(AN7,UO4)
Non_WC_distance_pymol(AN7,UN3)

#AN3
Non_WC_distance_pymol(AN3,GO6)
Non_WC_distance_pymol(AN3,GN1)
Non_WC_distance_pymol(AN3,GN2)
Non_WC_distance_pymol(AN3,GN3)
Non_WC_distance_pymol(AN3,GN7)

Non_WC_distance_pymol(AN3,CN3)
Non_WC_distance_pymol(AN3,CO2)

Non_WC_distance_pymol(AN3,AN6)
Non_WC_distance_pymol(AN3,AN7)
Non_WC_distance_pymol(AN3,AN9)

Non_WC_distance_pymol(AN3,UO4)
Non_WC_distance_pymol(AN3,UN3)


#UO4
Non_WC_distance_pymol(UO4,GN1)
Non_WC_distance_pymol(UO4,GN2)
Non_WC_distance_pymol(UO4,GN3)
Non_WC_distance_pymol(UO4,GN7)

Non_WC_distance_pymol(UO4,CN3)
Non_WC_distance_pymol(UO4,CO2)

Non_WC_distance_pymol(UO4,AN1)
Non_WC_distance_pymol(UO4,AN3)
Non_WC_distance_pymol(UO4,AN7)

#UN3
Non_WC_distance_pymol(UN3,GO6)
Non_WC_distance_pymol(UN3,GN1)
Non_WC_distance_pymol(UN3,GN2)
Non_WC_distance_pymol(UN3,GN3)
Non_WC_distance_pymol(UN3,GN7)

Non_WC_distance_pymol(UN3,CN3)
Non_WC_distance_pymol(UN3,CO2)

Non_WC_distance_pymol(UN3,AN6)
Non_WC_distance_pymol(UN3,AN3)
Non_WC_distance_pymol(UN3,AN7)

Non_WC_distance_pymol(UN3,UO4)
Non_WC_distance_pymol(UO4,UO2)

#CN4
Non_WC_distance_pymol(CN4,GN3)
Non_WC_distance_pymol(CN4,GN7)

Non_WC_distance_pymol(CN4,AN1)
Non_WC_distance_pymol(CN4,AN3)
Non_WC_distance_pymol(CN4,AN7)


Non_WC_distance_pymol(CN4,UO4)
Non_WC_distance_pymol(CN4,UN3)

#CN3
Non_WC_distance_pymol(CN3,GN3)
Non_WC_distance_pymol(CN3,GN7)

Non_WC_distance_pymol(CN3,AN6)
Non_WC_distance_pymol(CN3,AN1)
Non_WC_distance_pymol(CN3,AN3)
Non_WC_distance_pymol(CN3,AN7)

Non_WC_distance_pymol(CN3,UO4)
Non_WC_distance_pymol(CN3,UN3)

#CO2
Non_WC_distance_pymol(CO2,GO6)
Non_WC_distance_pymol(CO2,GN3)
Non_WC_distance_pymol(CO2,GN7)

Non_WC_distance_pymol(CO2,AN6)
Non_WC_distance_pymol(CO2,AN1)
Non_WC_distance_pymol(CO2,AN3)
Non_WC_distance_pymol(CO2,AN7)


Non_WC_distance_pymol(CO2,UO4)
Non_WC_distance_pymol(CO2,UN3)

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


#Done

'''

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('PyMIV', run_plugin_gui)

# To create global reference of the dialog variables
dialog = None

# to give filename of the UI file
uifile = os.path.join(os.path.dirname(__file__), 'pymolGUI.ui')

# To load the UI dialog
form = loadUi(uifile, dialog)


def run_plugin_gui():
    '''
    Open the custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()


def make_dialog():
    '''
    Make the dialog window
    '''
    # To create a new UI window
    dialog = QtWidgets.QDialog()

    # To populate the Window from our .ui file
    uifile = os.path.join(os.path.dirname(__file__), 'pymolGUI.ui')
    form = loadUi(uifile, dialog)
   

    def run():
        '''
        Run the actions in the dialog window
        '''
        # retreive PDB file data
        pdb_file = form.lineEdit.text()
        calc_disulfide(pdb_file)

        # To debug code
        print('User Entered Filename:', pdb_file)

    # To connect clicking buttons to a value, text or command
    form.disulfideFinder.clicked.connect(run)
    form.calculateMW.clicked.connect(run)
    form.wcAndNonWC.clicked.connect(run)
    form.hydrogenBond.clicked.connect(run)

    return dialog