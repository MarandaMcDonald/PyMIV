#Driver code
import re
from calc_interaction import *

calc_disulfide("PDB_Files/1fdl.pdb")
#alpha_helice("PDB_Files/1kx5.pdb")

calc_wc_nwc("PDB_Files/1kx5.pdb")
#calc_wc_nwc("PDB_Files/1z43.pdb")

#calc_peptide_mw("PDB_Files/1fdl.pdb")

#fastaSeqList=[]
#output_fasta("PDB_Files/1fdl.pdb", fastaSeqList)

#end_to_end_dist("PDB_Files/1fdl.pdb")
