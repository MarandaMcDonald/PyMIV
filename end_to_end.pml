load PDB_Files/1fdl.pdb
remove resn hoh
dist end_to_end, /1fdl//L/P`1/CA ,/1fdl//Y/U`129/CA
show sticks, /1fdl//L/P`1/CA 
show sticks, /1fdl//Y/U`129/CA
set dash_color, yellow, end_to_end
set dash_length, 0.2500
set dash_gap, 0.4
set dash_radius, .15
