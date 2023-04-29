#This is a series of functions for mathematical calculations to read into PDB files and output a .pml text to import into the PyMOL


import math

############################################################
#######################Core Functions#######################
############################################################

def PDB_file_read():
    '''
        This function will read into a PDB file format using the readlines function

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

def atom_finder(base, atom, mylist):
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

def residue_finder(pdblist,resn="CYS", atom="S", mylist="CYSlistUnsorted"):
    '''
    This function will read PDB file and search for a given amino acid resiue, and the atomic element specified. This atom will be appended to a list given

    **Parameters**

    PDBfile: *list*
        A list of strings in PDB file format
    **Returns**

        Sorted list by acsending resiue number position
    '''
    for line in pdblist:
        if (line[0:4]=="ATOM") and (line[17:20]==resn) and (line[77:78]==atom):
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

############################################################
#####################Detect Sulfide Bonds###################
############################################################

def calc_disulfide():
    '''
        This function will calculate any disulfide bonds in a PDB file and display the 

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
    ##function to find CYS
    CYSlistUnsorted=[]
    
    residue_finder(pdblist,"CYS", CYSlistUnsorted)

    CYSlist=sorted_pdb(CYSlistUnsorted)
    #To print the number of CYS residues

    print()
    print("There are",len(CYSlist), "CYS residues")
    #test extractxyz works and print the xyz coordinates for each CYS atoms
    #for line in CYSlist:
    #    temp=extractxyz(line)
    #    print("CYS Res{}, x={}, y={}, z={}".format(line[22:29],temp[0],temp[1],temp[2]))


    CYSdistanceslist=[]

    cyskeys= []
    #values are the bond distance
    cysvalues = []
    cysdict={}
    ### here loop through the PDB file, and if residues are sequential, then calc distances and
    # store in a new list
    for i in range(len(CYSlist)):
        for z in range(i):
            CYSdistanceslist.append(calcdistance(extractxyz(CYSlist[i]), extractxyz(CYSlist[z])))

    #to detect if a potential disulfide bond is within the accpetable 2.00 angstroms plus or minus 0.05
    def find_s_bond(mylist):
        for i in mylist:
            if 2.05 >i > 1.95:
                print("{}, {}".format(i[17:29],i))

    #for i in CYSdistanceslist:
    #    print(i)



    StoSlist=[]
    def s_to_s(mylist):
        for i in range(len(mylist)):
            for z in range(i):
                StoSlist.append((mylist[i]) + (mylist[z]))
    #print("Here are the potential S-S bond interactions")
    s_to_s(CYSlist)
    #keys are the Cys to Cys residue numbers

    #for line in StoSlist:
    #            print(line[17:27], "to", line[97:110], "S-S Distance:")

    #make list of cys distances
    for i in CYSdistanceslist:
        cysvalues.append(i)

    #make list of cys to cys combinations
    for i in StoSlist:
        cyskeys.append(i)
    #to make a dictionary matching the residue combination to the bond distance
    for i in range(len(cyskeys)):
        cysdict[cyskeys[i]] = cysvalues[i]


    #to print out the CYS RESN to CYS RESN combination
    #for i in cysdict.keys():
    #    print(i[17:27], "to", i[97:110], "S-S Distance:")

    cysValueToKey = {i for i in cysdict if 2.05 >cysdict[i] > 1.95}
    trueCysBondlist = list(cysValueToKey)
    trueDistancelist = list()
    for i in cysdict.values():
            if 2.05 > i> 1.95:
            trueDistancelist.append(i)

    def true_cys(trueCysBondList):
        for i in trueCysBondlist:
            print(i[17:27], "to", i[97:110], "S-S Distance:",)
    def true_distance(trueDistancelist):
        for i in trueDistancelist:
            print(i)

    #print(("{}, {}".format(true_cys(trueCysBondlist), true_distance(trueDistancelist))))

    joinedList = "\n".join("{} {}".format(x, y) for x, y in zip(trueCysBondlist, "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"))

    print("There are", len(trueDistancelist), "disulfide bonds")
    print()
    print("DISULFDE BONDS ( 2 ± 0.05 Å )")

    #print(joinedList)

    for i in trueCysBondlist:
        print(i[17:27], "---", i[97:107],)
    print()
    print("Thanks for using me!")

############################################################
###############WC and Non-WC Nucleic Acid Interactions######
############################################################
def 
'''
    This function will encrypt a message that is passed to it

    **Parameters**

    message: *str*
        The given plain text message that will be encrypted

    N: *int*
        The value for N

    E: *int*
        The value for E

    **Returns**
    
        None
    '''
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
(atom_finder("G","O6 ",GO6))
(atom_finder("G","N1 ",GN1))
(atom_finder("G","N2 ",GN2))

(atom_finder("C","N4 ",CN4))
(atom_finder("C","N3 ",CN3))
(atom_finder("C","O2 ",CO2))

(atom_finder("A","N6 ",AN6))
(atom_finder("A","N1 ",AN1))

(atom_finder("U","O4 ",UO4))
(atom_finder("U","N3 ",UN3))   

#Non-WC atoms
(atom_finder("G","N3 ",GN3))
(atom_finder("G","N9 ",GN9))
(atom_finder("G","N7 ",GN7))


(atom_finder("A","N7 ",AN7))
(atom_finder("A","N9 ",AN9))
(atom_finder("A","N3 ",AN3))

(atom_finder("U","O2 ",UO2))


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
def bond_maker(list1, list2):
    for i in range(len(list1)):
        for z in range(len(list2)):
            if list1==GO6 and list2==CN4:
                print("done")
            else:
                Non_WC_distance_pymol(list1[i], list2[z])'''


def non_wc(list1, list2):
    for i in range(len(atomList)):
        for z in range(len(atomList)):
            if atomList[i]==GO6 and atomList[z]==CN4:
                print("skip")
            else:
                Non_WC_distance_pymol(atomList[i], atomList[z])

'''
bond_maker(Ghoog,Ghoog)
bond_maker(Ghoog,Ahoog)
bond_maker(Ghoog,Choog)
bond_maker(Ghoog,Uhoog)

bond_maker(Gwc,Uwc)
bond_maker(Gwc,Awc)

bond_maker(Cwc,Uwc)
bond_maker(Cwc,Uwc)'''

'''
Non_WC_distance_pymol(,GO6)
Non_WC_distance_pymol(,GN1)
Non_WC_distance_pymol(,GN2)
Non_WC_distance_pymol(,GN3)
Non_WC_distance_pymol(,GN9)
Non_WC_distance_pymol(,GN7)

Non_WC_distance_pymol(,CN3)
Non_WC_distance_pymol(,CO2)


Non_WC_distance_pymol(,AN6)
Non_WC_distance_pymol(,AN1)
Non_WC_distance_pymol(,AN3)
Non_WC_distance_pymol(,AN7)
Non_WC_distance_pymol(,AN9)

Non_WC_distance_pymol(,UO4)
Non_WC_distance_pymol(,UN3)
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