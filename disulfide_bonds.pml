load PDB_Files/1fdl.pdb
remove resn hoh
dist disulfide_bond, /1fdl//H/CYS`218/SG ,/1fdl//L/CYS`214/SG
show sticks, /1fdl//H/CYS`218/SG 
show sticks, /1fdl//L/CYS`214/SG
color atomic, /1fdl//H/CYS`218/SG 
color atomic, /1fdl//L/CYS`214/SG
dist disulfide_bond, /1fdl//L/CYS`88/SG ,/1fdl//L/CYS`23/SG
show sticks, /1fdl//L/CYS`88/SG 
show sticks, /1fdl//L/CYS`23/SG
color atomic, /1fdl//L/CYS`88/SG 
color atomic, /1fdl//L/CYS`23/SG
dist disulfide_bond, /1fdl//L/CYS`194/SG ,/1fdl//L/CYS`134/SG
show sticks, /1fdl//L/CYS`194/SG 
show sticks, /1fdl//L/CYS`134/SG
color atomic, /1fdl//L/CYS`194/SG 
color atomic, /1fdl//L/CYS`134/SG
dist disulfide_bond, /1fdl//Y/CYS`94/SG ,/1fdl//Y/CYS`76/SG
show sticks, /1fdl//Y/CYS`94/SG 
show sticks, /1fdl//Y/CYS`76/SG
color atomic, /1fdl//Y/CYS`94/SG 
color atomic, /1fdl//Y/CYS`76/SG
dist disulfide_bond, /1fdl//Y/CYS`115/SG ,/1fdl//Y/CYS`30/SG
show sticks, /1fdl//Y/CYS`115/SG 
show sticks, /1fdl//Y/CYS`30/SG
color atomic, /1fdl//Y/CYS`115/SG 
color atomic, /1fdl//Y/CYS`30/SG
dist disulfide_bond, /1fdl//Y/CYS`127/SG ,/1fdl//Y/CYS`6/SG
show sticks, /1fdl//Y/CYS`127/SG 
show sticks, /1fdl//Y/CYS`6/SG
color atomic, /1fdl//Y/CYS`127/SG 
color atomic, /1fdl//Y/CYS`6/SG
dist disulfide_bond, /1fdl//H/CYS`95/SG ,/1fdl//H/CYS`22/SG
show sticks, /1fdl//H/CYS`95/SG 
show sticks, /1fdl//H/CYS`22/SG
color atomic, /1fdl//H/CYS`95/SG 
color atomic, /1fdl//H/CYS`22/SG
dist disulfide_bond, /1fdl//Y/CYS`80/SG ,/1fdl//Y/CYS`64/SG
show sticks, /1fdl//Y/CYS`80/SG 
show sticks, /1fdl//Y/CYS`64/SG
color atomic, /1fdl//Y/CYS`80/SG 
color atomic, /1fdl//Y/CYS`64/SG
dist disulfide_bond, /1fdl//H/CYS`198/SG ,/1fdl//H/CYS`143/SG
show sticks, /1fdl//H/CYS`198/SG 
show sticks, /1fdl//H/CYS`143/SG
color atomic, /1fdl//H/CYS`198/SG 
color atomic, /1fdl//H/CYS`143/SG

hide labels, disulfide_bond
set dash_length, 0.2500
set dash_gap, 0.4
set dash_radius, .15
