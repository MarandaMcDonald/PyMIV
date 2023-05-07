load PDB_Files/1fdl.pdb
remove resn hoh
dist disulfide_bond, /1fdl//Y/CYS`94/SG ,/1fdl//Y/CYS`76/SG
show sticks, /1fdl//Y/CYS`94
show sticks, /1fdl//Y/CYS`76
color atomic, /1fdl//Y/CYS`94
color atomic, /1fdl//Y/CYS`76
dist disulfide_bond, /1fdl//H/CYS`95/SG ,/1fdl//H/CYS`22/SG
show sticks, /1fdl//H/CYS`95
show sticks, /1fdl//H/CYS`22
color atomic, /1fdl//H/CYS`95
color atomic, /1fdl//H/CYS`22
dist disulfide_bond, /1fdl//Y/CYS`127/SG ,/1fdl//Y/CYS`6/SG
show sticks, /1fdl//Y/CYS`127
show sticks, /1fdl//Y/CYS`6
color atomic, /1fdl//Y/CYS`127
color atomic, /1fdl//Y/CYS`6
dist disulfide_bond, /1fdl//L/CYS`194/SG ,/1fdl//L/CYS`134/SG
show sticks, /1fdl//L/CYS`194
show sticks, /1fdl//L/CYS`134
color atomic, /1fdl//L/CYS`194
color atomic, /1fdl//L/CYS`134
dist disulfide_bond, /1fdl//Y/CYS`80/SG ,/1fdl//Y/CYS`64/SG
show sticks, /1fdl//Y/CYS`80
show sticks, /1fdl//Y/CYS`64
color atomic, /1fdl//Y/CYS`80
color atomic, /1fdl//Y/CYS`64
dist disulfide_bond, /1fdl//H/CYS`198/SG ,/1fdl//H/CYS`143/SG
show sticks, /1fdl//H/CYS`198
show sticks, /1fdl//H/CYS`143
color atomic, /1fdl//H/CYS`198
color atomic, /1fdl//H/CYS`143
dist disulfide_bond, /1fdl//H/CYS`218/SG ,/1fdl//L/CYS`214/SG
show sticks, /1fdl//H/CYS`218
show sticks, /1fdl//L/CYS`214
color atomic, /1fdl//H/CYS`218
color atomic, /1fdl//L/CYS`214
dist disulfide_bond, /1fdl//Y/CYS`115/SG ,/1fdl//Y/CYS`30/SG
show sticks, /1fdl//Y/CYS`115
show sticks, /1fdl//Y/CYS`30
color atomic, /1fdl//Y/CYS`115
color atomic, /1fdl//Y/CYS`30
dist disulfide_bond, /1fdl//L/CYS`88/SG ,/1fdl//L/CYS`23/SG
show sticks, /1fdl//L/CYS`88
show sticks, /1fdl//L/CYS`23
color atomic, /1fdl//L/CYS`88
color atomic, /1fdl//L/CYS`23

hide labels, disulfide_bond
set dash_length, 0.2500
set dash_gap, 0.4
set dash_radius, .15
