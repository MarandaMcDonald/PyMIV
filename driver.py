#Driver code

from calc_interaction import *

calc_disulfide("PDB_Files/1fdl.pdb")
alpha_helice("PDB_Files/1njg.pdb")

calc_wc_nwc("PDB_Files/1z43.pdb")

calc_peptide_mw("PDB_Files/1fdl.pdb")

fastaSeqList=[]
output_fasta("PDB_Files/1fdl.pdb", fastaSeqList)

end_to_end_dist("PDB_Files/1fdl.pdb")

test="(['/Users/marandamacpro/Downloads/PyMIV/PDB_Files/1bhm.pdb'], 'All Files (*)')"

#clean_path=clean_file_path_one(clean_file_path_two)()
#print(clean_path(test))

print(test[3:-20])