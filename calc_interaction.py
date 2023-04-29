#This is a series of functions for mathematical calculations to read into PDB files and output a .pml text to import into the PyMOL



#run "get_bonds.pml" in pymol with the 1z43.pdb structure
import math

#to remove empty spaces in strings
def remove(string):
    return string.replace(" ", "")
    
def extractxyz(PDBline):
  ### from expected PDB line, return x y z as a list of floats
  return([float(PDBline[29:38]),float(PDBline[38:46]),float(PDBline[46:54])])

def calcdistance(list1,list2):
  #expects two lists of floats, x1, y1, z1 and x2, y2, z2
  return(math.sqrt((list1[0]-list2[0])**2 + (list1[1]-list2[1])**2 + (list1[2]-list2[2])**2))
list1=[]
def atom_finder(base, atom, mylist):
    for line in pdblist:
        if (line[0:4]=="ATOM"):
            if (line[19:20]==base):
                if (line[13:16]==atom):
                    mylist.append(line)



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
GO6=[]
GN1=[]
GN2=[]

CN4=[]
CN3=[]
CO2=[]

AN6=[]
AN1=[]

UO4=[]
UN3=[]

#lists of Non-WC bonding atoms
GN3=[]
GN9=[]
GN7=[]


AN7=[]
AN9=[]
AN3=[]

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